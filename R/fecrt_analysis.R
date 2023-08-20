#' An example R6 class
#'
#' @param value an initialisation value
#' @param new_value a new value to use
#'
#' @importFrom R6 R6Class
#' @importFrom tools file_ext
#' @importFrom readxl read_excel excel_sheets
#' @importFrom readr read_csv read_csv2
#' @importFrom stringr str_c
#'
#' @examples
#' example_instance <- FecrtAnalysis$new(shiny = TRUE)
#' example_instance
#'
#' @export
FecrtAnalysis <- R6Class("FecrtAnalysis",

	public = list(

		initialize = function(shiny = FALSE){
			private$shiny <- shiny
		},

		reset = function(){
			private$current_type <- ""
			private$data <- list()
		},

		add_data_files = function(files, names=files){

			seq_len(length(files)) |>
				as.list() |>
				lapply(function(i){
					ext <- file_ext(names)
					ext <- ifelse(length(ext)==0L, "", ext)

					if(ext %in% "csv"){
						rv <- self$add_data_csv(files[i], name=names[i])
					}else if(ext %in% c("xlsx","xls")){
						rv <- self$add_data_excel(files[i], name=names[i])
					}else{
						rv <- str_c(names[i], ": ERROR : Invalid file (please upload only .csv, .xls or .xlsx files)")
					}
					rv
				}) |>
				unlist() ->
				messages

			messages
		},

		add_data_excel = function(file, name=file){
			stopifnot(length(file)==1)
			stopifnot(file.exists(file))

			## Add data from each sheet:
			excel_sheets(file) |>
				as.list() |>
				lapply(function(s){
					read_excel(file, s, .name_repair="unique") |>
						self$import_data(name=s) |>
						c(str_c(name, " - ")) |>
						rev() |>
						str_c(collapse="")
				}) |>
				unlist()
		},

		add_data_csv = function(file, name=file){
			stopifnot(length(file)==1)
			stopifnot(file.exists(file))

			## Try to auto-read csv or csv2:
			data <- readr::read_csv(file, show_col_types = FALSE)
			if(ncol(data)==1L){
				data <- readr::read_csv2(file, show_col_types = FALSE)
			}

			self$import_data(data, name)
		},

		import_data = function(data, name){
			stopifnot(is.data.frame(data))
			stopifnot(is.character(name), length(name)==1L, !is.na(name), name!="")
			if(name %in% sapply(private$data, function(x) x$name)){
				return(str_c(name, ": ERROR : dataset with matching name previously imported"))
			}

			if(nrow(data)==0L){
				return(str_c(name, ": ERROR : zero rows of data"))
			}

			names(data) |>
				str_replace_all(" ", "") |>
				str_replace_all("-", "") |>
				tolower() ->
				cols

			types <- logical(3)
			names(types) <- c("paired_wide","unpaired_wide","unpaired_long")

			if(all(c("pretreatment","posttreatment") %in% cols)){
				types["paired_wide"] <- TRUE
			}
			if(all(c("treatment","control") %in% cols)){
				types["unpaired_wide"] <- TRUE
			}
			if(all(c("group","epg") %in% cols)){
				types["unpaired_long"] <- TRUE
			}

			if(sum(types)==0L){
				return(str_c(name, ": ERROR : column names do not match accepted data formats"))
			}else if(sum(types)>1L){
				return(str_c(name, ": ERROR : multiple column names detected in different formats"))
			}

			if(types["paired_wide"]){
				if(private$current_type=="") private$current_type <- "paired"
				if(private$current_type!="paired"){
					return(str_c(name, ": ERROR : data format does not match that of previously imported datasets"))
				}

				m1 <- which(cols %in% "pretreatment")
				m2 <- which(cols %in% "posttreatment")
				if(length(m1)>1 || length(m2)>1){
					return(str_c(name, ": ERROR : multiple columns matching PreTreatment and/or PostTreatment"))
				}
				d1 <- data[[m1]]
				d2 <- data[[m2]]

				if(all(is.na(d1))){
					return(str_c(name, ": ERROR : values in the PreTreatment column are all missing"))
				}
				if(!is.numeric(d1)){
					return(str_c(name, ": ERROR : some values in the PreTreatment column are not numbers"))
				}
				if(any(d1 < 0, na.rm=TRUE)){
					return(str_c(name, ": ERROR : one or more values in the PreTreatment column are negative"))
				}
				if(all(is.na(d2))){
					return(str_c(name, ": ERROR : values in the PreTreatment column are all missing"))
				}
				if(!is.numeric(d2)){
					return(str_c(name, ": ERROR : some values in the PostTreatment column are not numbers"))
				}
				if(any(d2 < 0, na.rm=TRUE)){
					return(str_c(name, ": ERROR : one or more values in the PreTreatment column are negative"))
				}

				newdata <- list()
				private$data <- c(private$data, list(list(name=name, type="paired", PreTreatment=d1, PostTreatment=d2)))
				sumobs <- sum(c(!is.na(d1), !is.na(d2)))

			}else if(types["unpaired_wide"]){
				if(private$current_type=="") private$current_type <- "unpaired"
				if(private$current_type!="unpaired"){
					return(str_c(name, ": ERROR : data format does not match that of previously imported datasets"))
				}

				m1 <- which(cols %in% "control")
				m2 <- which(cols %in% "treatment")
				if(length(m1)>1 || length(m2)>1){
					return(str_c(name, ": ERROR : multiple columns matching Control and/or Treatment"))
				}
				d1 <- data[[m1]]
				d2 <- data[[m2]]

				if(all(is.na(d1))){
					return(str_c(name, ": ERROR : values in the Control column are all missing"))
				}
				if(!is.numeric(d1)){
					return(str_c(name, ": ERROR : some values in the Control column are not numbers"))
				}
				if(any(d1 < 0, na.rm=TRUE)){
					return(str_c(name, ": ERROR : one or more values in the Control column are negative"))
				}
				if(all(is.na(d2))){
					return(str_c(name, ": ERROR : values in the Treatment column are all missing"))
				}
				if(!is.numeric(d2)){
					return(str_c(name, ": ERROR : some values in the Treatment column are not numbers"))
				}
				if(any(d2 < 0, na.rm=TRUE)){
					return(str_c(name, ": ERROR : one or more values in the Treatment column are negative"))
				}

				newdata <- list()
				private$data <- c(private$data, list(list(name=name, type="unpaired", Control=d1, Treatment=d2)))
				sumobs <- sum(c(!is.na(d1), !is.na(d2)))

			}else if(types["unpaired_long"]){
				if(private$current_type=="") private$current_type <- "unpaired"
				if(private$current_type!="unpaired"){
					return(str_c(name, ": ERROR : data format does not match that of previously imported datasets"))
				}

				m1 <- which(cols %in% "group")
				m2 <- which(cols %in% "epg")
				if(length(m1)>1 || length(m2)>1){
					return(str_c(name, ": ERROR : multiple columns matching Group and/or EPG"))
				}
				d1 <- data[[m1]]
				d2 <- data[[m2]]

				if(any(is.na(d1))){
					return(str_c(name, ": ERROR : some values in the Group column are missing"))
				}
				if(!all(tolower(d1) %in% c("control","treatment"))){
					return(str_c(name, ": ERROR : some values in the Group column do not match 'Control' or 'Treatment'"))
				}

				if(all(is.na(d2))){
					return(str_c(name, ": ERROR : values in the EPG column are all missing"))
				}
				if(!is.numeric(d2)){
					return(str_c(name, ": ERROR : some values in the EPG column are not numbers"))
				}
				if(any(d2 < 0, na.rm=TRUE)){
					return(str_c(name, ": ERROR : one or more values in the EPG column are negative"))
				}

				ctl <- as.numeric(na.omit(d2[tolower(d1)=="control"]))
				txt <- as.numeric(na.omit(d2[tolower(d1)=="treatment"]))

				if(length(ctl)==0L){
					return(str_c(name, ": ERROR : zero Control animals with non-missing EPG"))
				}
				if(length(txt)==0L){
					return(str_c(name, ": ERROR : zero Treatment animals with non-missing EPG"))
				}

				newdata <- list()
				private$data <- c(private$data, list(list(name=name, type="unpaired", Control=ctl, Treatment=txt)))
				sumobs <- length(ctl)+length(txt)

			}else{
				stop("Unrecognised type")
			}

			return(str_c(name, ": OK : dataset imported (N=", sumobs, ")"))
		},

		get_data = function(){
			private$data
		},

		set_parameters_guidelines = function(guideline, version, mf_1, mf_2=mf_1, region=NA_character_, identifier=NA_character_){
			private$parameters <- list(guideline=guideline, version=version, mf_1=mf_1, mf_2=mf_2, region=region, identifier=identifier)
		},

		set_parameters = function(){
			## TODO: custom target efficacies etc
		},

		run_analysis = function(){
			stopifnot(!identical(private$data, list()), !identical(private$parameters, list()))

			stopifnot(private$current_type%in%c("paired","unpaired"))
			paired <- private$current_type=="paired"

			stopifnot(private$parmeters$version%in%c("clinical","research"))
			version <- private$parameters$version

			TE <- switch(private$parameters$guideline,
				"cattle" = 0.99,
				"sheep" = 0.99,
				"goats" = 0.99,
				"cyath_ml" = 0.999,
				"cyath_bz" = 0.99,
				"cyath_pyr" = 0.98,
				"foals" = 0.999,
				"pig_bz" = 0.99,
				"pig_ivm" = 0.95,
				NA_real_
			)

			if(version=="research"){
				TL <- switch(private$parameters$guideline,
					"cattle" = 0.95,
					"sheep" = 0.95,
					"goats" = 0.95,
					"cyath_ml" = 0.96,
					"cyath_bz" = 0.95,
					"cyath_pyr" = 0.88,
					"foals" = 0.95,
					"pig_bz" = 0.93,
					"pig_ivm" = 0.85,
					NA_real_
				)
			}else{
				TL <- switch(private$parameters$guideline,
					"cattle" = 0.90,
					"sheep" = 0.90,
					"goats" = 0.90,
					"cyath_ml" = 0.92,
					"cyath_bz" = 0.90,
					"cyath_pyr" = 0.80,
					"foals" = 0.90,
					"pig_bz" = 0.88,
					"pig_ivm" = 0.80,
					NA_real_
				)
			}
			stopifnot(!is.na(TE), !is.na(TL))

			k1 <- switch(private$parameters$guideline,
				"cattle" = 1.1,
				"sheep" = 1.4,
				"goats" = 3.8,
				"cyath_ml" = 1.3,
				"cyath_bz" = 1.3,
				"cyath_pyr" = 1.3,
				"foals" = 1.1,
				"pig_bz" = 0.8,
				"pig_ivm" = 0.8,
				NA_real_
			)
			k2 <- switch(private$parameters$guideline,
				"cattle" = 0.4,
				"sheep" = 0.8,
				"goats" = 2.4,
				"cyath_ml" = 0.5,
				"cyath_bz" = 0.5,
				"cyath_pyr" = 0.5,
				"foals" = 1.0,
				"pig_bz" = 1.2,
				"pig_ivm" = 1.2,
				NA_real_
			)
			kc <- switch(private$parameters$guideline,
				"cattle" = 0.1,
				"sheep" = 0.5,
				"goats" = 0.5,
				"cyath_ml" = 0.5,
				"cyath_bz" = 0.5,
				"cyath_pyr" = 0.5,
				"foals" = 0.2,
				"pig_bz" = 0.4,
				"pig_ivm" = 0.4,
				NA_real_
			)
			stopifnot(!is.na(k1), !is.na(k2), !is.na(kc))
			if(paired){
				k1 <- k1 / (1-kc)
				k2 <- k2 / (1-kc)
			}

			## Scale by mf for BNB and BNB_KnownKs:
			mf_1 <- private$parameters$mf_1
			mf_2 <- private$parameters$mf_2
			mf_r <- mf_1/mf_2
			TE_r <- 1- ((1-TE)*mf_r)
			TL_r <- 1- ((1-TL)*mf_r)

			for(i in self$n_data){
				if(paired){
					d1 <- private$data[[i]]$PreTreatment
					d2 <- private$data[[i]]$PostTreatment
				}else{
					d1 <- private$data[[i]]$PreTreatment
					d2 <- private$data[[i]]$PostTreatment
				}

				## Results for Gamma and WAAVP ignoring mf:
				d2 <- rnbinom(20,1,mu=1)*50
				results_nonpar <- efficacy_analysis(d1, d2, paired=paired, T_I=TE, T_A=TL, alpha=0.05)

				## Results for BNB and BNB_KnownKs:
				results_bnbs <- efficacy_analysis(d1/mf_1, d2/mf_2, paired=paired, T_I=TE_r, T_A=TL_r, alpha=0.05, known_ks = c(k1,k2))

				results <- bind_rows(
					results_nonpar |> filter(Method %in% c("Gamma","WAAVP")),
					results_bnbs |> filter(Method %in% c("BNB","BNB_KnownKs"))
				)
				results

				## NB: elements 4 and 5 are the real ks, 1 and 2 are adjusted!
				estks <- estimate_k(d1/mf_1, d2/mf_2, paired=paired)

				## If N<5 then prefer BNB_Knownks
				## If N>=5 and 3 or more non-zero counts (pre- and post-) then prefer Gamma
				## If N>=5 and <3 non-zero counts then prefer BNB
				## If any count %% 1 != 0 then don't show BNB or BNBKnownKs
				## Show WAABP as well


				## Code for other one:
				res <- fecrt_sample_size(EPGrange=50*20, EDT=50, mulength=1, margin=0.04, iters=10000L)
				# EPGrange=epgrange, mulength=length(epgrange), EDT=1, target=using$Target[1]/100, margin=unique(using$Delta/100), k1=using$k1[1], k2=using$k2[1], cor=using$cor[1], iters=iters
				ggplot(res$smat, aes(x=N, y=power, col=Hypothesis)) +
					geom_line() +
					geom_point() +
					scale_x_continuous(trans="log10")

			}

		}

	),

	active = list(
		n_data = function() length(private$data),
		design = function() private$current_type
	),

	private = list(
		shiny = FALSE,
		current_type = "",
		data = list(),
		parameters = list()
	)

)

if(FALSE){
	test <- FecrtAnalysis$new(shiny = TRUE)
	test$add_data_files("~/Downloads/file3.xlsx", "file3.xlsx")
	test$n_data
	test$import_data(data.frame(PreTreatment=1:10, PostTreatment=1:10), "test")
	test$import_data(data.frame(PreTreatment=1:10, PostTreatment=1:10), "test2")
	test$import_data(data.frame(Control=1:10, Treatment=1:10), "test3")
	test$import_data(data.frame(Group=c("Treatment","Control"), EPG=1:2), "test4")
	test$n_data

	ll <- letters
	ll[1] <- "a.csv"
	ll[2] <- "b.xlsx"
	rv <- test$add_data_files(ll)
	length(rv)
	head(rv, 3)
}
