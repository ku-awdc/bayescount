tempdir <- tempfile(pattern = "fecrt_analysis", tmpdir = tempdir(check=TRUE))
dir.create(tempdir)
on.exit(unlink(tempdir, recursive=TRUE))

fecrt_analysis <- FecrtAnalysis$new(shiny = TRUE)

function(input, output, session) {

  rv <- reactiveValues(
    data_paired = data.frame(PreTreatment=integer(0),PostTreatment=integer(0)),
    data_txt = data.frame(),
    data_ctl = data.frame(),
    data_demo = data.frame(),
    status = 0L,
    status_feedback = status_feedback[1],
    result_summary = ""
  )
  ## Note: statuses are:
  ## 0: nothing
  ## 1: N input (direct input only)
  ## 2: data input or uploaded
  ## 3: parameters selected and verified (multiplication factor != 0 etc)
  ## 4: summary results run
  ## 5: full results and markdown files generated

  settings <- reactiveValues(
    entryType = "upload",
    design = "paired",
    directN_paired = 20,
    directN_txt = 20,
    directN_ctl = 20
  )

  observe({
    req(input$entryType)

    if(input$entryType == "demo"){
      updateSelectInput(session, "design", selected="paired")
      updateNumericInput(session, "directN_paired", value=20L)

      settings$entryType <- "demo"
      settings$design <- "paired"
      settings$directN_paired <- 20L
      settings$directN_txt <- input$directN_txt
      settings$directN_ctl <- input$directN_ctl

      rv$data_demo <- data.frame(
        PreTreatment = as.integer(rnbinom(20, 1, mu=20)*50),
        PostTreatment = as.integer(rnbinom(20, 1, mu=1)*50)
      )

      waavp_choices_demo <- waavp_choices
      waavp_choices_demo$Swine <- NULL
      waavp_choices_demo$Equine <- waavp_choices_demo$Equine[1:3]

      updateSelectInput(session, "parameterType", selected="waavp")
      updateSelectInput(session, "waavpSetup", selected=sample(unlist(waavp_choices_demo),1))
      updateNumericInput(session, "mf_pre", value=50)
      updateCheckboxInput(session, "mfp_fixed", value=TRUE)
      updateTextInput(session, "region", value=sample(
        c("Westeros","Narnia","Discworld","Middle-earth","Aiur","Azeroth"), 1
      ))
      updateTextInput(session, "identifier", value="Demonstration")

      print("DEMO")
      rv$status <- 3
    }
  })

  data_file <- reactive({
    req(input$dataFile)

    fecrt_analysis$reset()
    msg <- fecrt_analysis$add_data_files(input$dataFile[,"datapath"], names=input$dataFile[,"name"])

    if(fecrt_analysis$n_data==0L){
      rv$status <- 0
      rv$status_feedback <- status_feedback[1]
      validate(str_c(c("No valid datasets found - please correct the following errors and try again:", str_c("\t", msg)), collapse="\n"))
    }

    msg <- str_c(c(str_c("A total of ", fecrt_analysis$n_data, " datasets were uploaded successfully with the following notes:"), str_c("\t", msg)), collapse="<br>")

    updateSelectInput(session, "design", selected = fecrt_analysis$design)
    settings$design <- fecrt_analysis$design

    ## Note: this is needed as otherwise switching back to direct breaks design:
    updateRadioButtons(session, "entryType", choices = c(`File upload` = "file"))

    rv$status <- 2
    rv$status_feedback <- status_feedback[2]

    msg
  })

  output$upload_feedback <- renderText({
    data_file()
  })

  observe({
    input$design
    rv$design <- input$design
    settings$design <- input$design
    rv$status <- 0L
  })

  output$status <- renderText(rv$status)
  outputOptions(output, 'status', suspendWhenHidden=FALSE)

  fluidPage({
    output$status_feedback <- renderText(rv$status_feedback)
    output$result_summary <- renderText(rv$result_summary)
  })


  ## Initialise data to go from status 0 to 1:
  observeEvent(input$initialise_data, {

    settings$design <- input$design
    settings$directN_paired <- input$directN_paired
    settings$directN_ctl <- input$directN_ctl
    settings$directN_txt <- input$directN_txt
    settings$entryType <- input$entryType

    rv$data_paired <- data.frame(
      `PreTreatment` = repnull(input$directN_paired),
      `PostTreatment` = repnull(input$directN_paired)
    )
    rv$data_ctl <- data.frame(
      `Control` = repnull(input$directN_ctl)
    )
    rv$data_txt <- data.frame(
      `Treatment` = repnull(input$directN_txt)
    )
    rv$status <- 1L
  })

  ## Input data to go from status 1 to 2 (upload done elsewhere):
  observe({
    rv$status
    input$data_paired
    input$data_demo
    input$data_ctl
    input$directN_txt

    if(rv$status >= 1 && settings$entryType %in% c("direct","demo")){
      if(settings$design == "paired"){
        if(settings$entryType=="direct"){
          data_paired <- hot_to_r(input$data_paired)
        }else{
          data_paired <- hot_to_r(input$data_demo)
        }
        rv$status <- ifelse(
          !is.null(data_paired) &&
            sum(!is.na(data_paired[[1]])) > 0L &&
            sum(!is.na(data_paired[[2]])) > 0L,
          max(rv$status, 2), 1)
      }else if(settings$design == "unpaired"){
        stopifnot(settings$entryType=="direct")
        data_ctl <- hot_to_r(input$data_ctl)
        data_txt <- hot_to_r(input$data_txt)
        rv$status <- ifelse(
          !is.null(data_ctl) &&
            !is.null(data_txt) &&
            sum(!is.na(data_ctl)) > 0L &&
            sum(!is.na(data_txt)) > 0L,
          max(rv$status, 2), 1)
      }else{
        stop("Unrecognised settings$design")
      }
      if(rv$status==2){
        rv$status_feedback <- status_feedback[2]
      }else{
        rv$status_feedback <- status_feedback[1]
      }
    }
    ## Note: file input checking is done elsewhere!

    print(rv$status)
  })

  ## Input parameters to go from status 2 to 3:
  observe({
    rv$status
    input$parameterType
    input$waavpSetup
    input$mf_pre
    input$mf_post
    input$mfp_fixed
    input$mf_txt
    input$mf_ctl
    input$mfu_fixed

    if(rv$status >= 2){
      if(input$design=="paired"){
        if(is.na(input$mf_pre) || (!input$mfp_fixed && is.na(input$mf_post))){
          rv$status <- 2
        }else{
          rv$status <- max(rv$status, 3)
        }
      }else if(input$design=="unpaired"){
        if(is.na(input$mf_ctl) || (!input$mfu_fixed && is.na(input$mf_txt))){
          rv$status <- 2
        }else{
          rv$status <- max(rv$status, 3)
        }
      }else{
        stop("Unrecognised design")
      }
    }

    print(rv$status)
  })


  ## Settings that reset status to 0:
  observe({
    if(changed(c("entryType","design","directN_paired","directN_ctl","directN_txt"),input,settings)){
      ## If moving away from demo then reset everything:
      if(settings$entryType=="demo"){
        print("HARD RESET")
        updateSelectInput(session, "design", selected="paired")
        updateNumericInput(session, "directN_paired", value=20L)
        updateSelectInput(session, "parameterType", selected="waavp")
        updateSelectInput(session, "waavpSetup", selected=unlist(waavp_choices)[1])
        updateNumericInput(session, "mf_pre", value=NULL)
        updateCheckboxInput(session, "mfp_fixed", value=TRUE)
        updateTextInput(session, "region", value="")
        updateTextInput(session, "identifier", value="")
        settings$entryType <- input$entryType
      }

      rv$status <- 0L
    }
  })


  observe({
    output$data_paired <- renderRHandsontable({
      rhandsontable(rv$data_paired, rowHeaders=NULL, stretchH = "none")
    })
    output$data_demo <- renderRHandsontable({
      rhandsontable(rv$data_demo, rowHeaders=NULL, stretchH = "none")
    })
    output$data_txt <- renderRHandsontable({
      rhandsontable(rv$data_txt, rowHeaders=NULL, stretchH = "none")
    })
    output$data_ctl <- renderRHandsontable({
      rhandsontable(rv$data_ctl, rowHeaders=NULL, stretchH = "none")
    })
  })

  observe({

    parasite_choices <-
      if(input$hostSpecies == "INVALID"){
        "*SELECT HOST FIRST*"
      }else if(input$hostSpecies %in% c("cattle","sheep","goats")){
        c("*SELECT*", "nematodes", "other")
      }else if(input$hostSpecies %in% c("horse_adult","donk_adult")){
        c("*SELECT*", "strongyles", "other")
      }else{
        "ERROR - MISSING HOST"
      }

    anthelmintic_choices <-
      if(input$hostSpecies == "INVALID"){
        "*SELECT HOST FIRST*"
      }else if(input$hostSpecies %in% c("cattle","sheep","goats")){
        c("*SELECT*", "Benzimidazoles", "Other")
      }else if(input$hostSpecies %in% c("horse_adult","donk_adult")){
        c("*SELECT*", "Pyrantel", "Other")
      }else{
        "ERROR - MISSING HOST"
      }

    updateSelectInput(session,
      "parasiteSpecies",
      choices = parasite_choices,
      selected = ifelse(is.null(input$parasiteSpecies) || !input$parasiteSpecies %in% parasite_choices, parasite_choices[1], input$parasiteSpecies)
    )

    updateSelectInput(session,
      "anthelmintic",
      choices = anthelmintic_choices,
      selected = ifelse(is.null(input$anthelmintic) || !input$anthelmintic %in% anthelmintic_choices, anthelmintic_choices[1], input$anthelmintic)
    )

  })

  #output$host <- "other"
  output$footer <- renderText("<h4>Hello</h4>")
  #print(input$hostSpecies)
  #print(input)

  output$downloadReport <- downloadHandler(
    filename = function() {
      ext <- if(input$downloadType=="pdf"){
        "pdf"
      }else if(input$downloadType=="word"){
        "docx"
      }else{
        stop("Unrecognised download type")
      }
      paste(ifelse(identical(input$identifier,""), "report", input$identifier), ".", ext, sep = "")
    },
    content = function(file) {
      ext <- if(input$downloadType=="pdf"){
        "pdf"
      }else if(input$downloadType=="word"){
        "docx"
      }else{
        stop("Unrecognised download type")
      }
      success <- file.copy(file.path(tempdir, paste("report.",ext,sep="")),file)
      if(!success) stop("Error copying report")
    }
  )

  observeEvent(input$calculate, {
    withProgress(message = "Calculating...", value= 0, {

      ## Create R6 class object to check/store everything:

      for(i in 1:10){
        setProgress(i/10)
        Sys.sleep(0.1)
      }
    })

    rv$result_summary <- "Results will be displayed here...<br>Note: the report download links currently do not do anything!"
    rv$status <- 4
  })


}
