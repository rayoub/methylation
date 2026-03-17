
suppressWarnings(suppressPackageStartupMessages({
	library(shiny)
	library(DT)
	library(here)
	library(bslib)
}))

source(here::here("R", "files.R"))

# ********************************************************************************
# *** functions
# ********************************************************************************

getBatchIds <- function () {

	dir <- here::here("results", "lab")
	files <- list.files(dir, "*_betas.rds")
	batch_ids <- gsub("_betas.rds","", files)

	return (batch_ids)
}

# ********************************************************************************
# *** ui
# ********************************************************************************

ui <- fluidPage(
	titlePanel("Methylation CNS Tumor Classification"),
	tabsetPanel(
		tabPanel("Import data", 
			br(),
			selectInput("dataset", label = "Dataset", choices = ls("package:datasets")),
			verbatimTextOutput("summary"),
			tableOutput("table")
		),
		tabPanel("View data",
			br(),
			selectInput("selectBatch", label = "Select Batch", choices = getBatchIds()),
			DTOutput("sampleTable")
		)
	)
)

# ********************************************************************************
# *** server
# ********************************************************************************

server <- function (input, output, session) {

	# *** Import tab

	dataset <- reactive({
		get(input$dataset, "package:datasets")
	})
	
	output$summary <- renderPrint({
		summary(dataset())
	})

	output$table <- renderTable({
		dataset()
	})

	# *** View tab

	r_scores <- reactive({
		message("r_scores triggered")
		scores <-loadLabData(input$selectBatch, "scores")
		scores_df <- as.data.frame(scores, row.names = rownames(scores))
		scores_t <- tibble::rownames_to_column(scores_df, var = "sample_id")
	})

	output$sampleTable <- renderDT({

		scores_mc <- r_scores() |>
			pivot_longer(
				cols = !(sample_id),
				names_to = "mc",
				values_to = "mc_score"
			) |>
			dplyr::group_by(sample_id) |>
			dplyr::slice_max(mc_score) 

		scores_mcf <- r_scores() |>
			pivot_longer(
				cols = !(sample_id),
				names_to = "mc",
				values_to = "mc_score"
			) |>
			dplyr::mutate(
				mcf = mcf_lookup(mc),
				.after = mc
			) |>
			dplyr::group_by(sample_id, mcf) |>
			dplyr::summarize(
				mcf_score = sum(mc_score),
				.groups = "drop_last"
			) |>
			dplyr::slice_max(order_by = mcf_score)

		scores_j <- scores_mc |> 
			dplyr::inner_join(scores_mcf, join_by(sample_id)) |>
			dplyr::mutate(
				mc_descr = MC[mc]
			) 

		datatable(scores_j, 
			colnames = c("Sample ID", "Methylation Class", "MC Score", "Methylation Class Family","MCF Score", "Methylation Class Description"),
			rowname = FALSE,
			selection = list(mode = "single", selected = 1),
			options = list(
					lengthChange = FALSE,
					pageLength = 10
				)
			) |>
			formatPercentage(columns = c("mc_score", "mcf_score"), digits = 2) 

	}, server = FALSE)


#	observeEvent(r_scores(), {
#		message("event triggered")
#		freezeReactiveValue(input, "selectSample")
#		updateSelectInput(inputId = "selectSample", choices = r_scores()$sample_id)
#	})

#	r_scores_max <- reactive({
#		message("r_scores_max trigged")
#		r_scores() |>
#			pivot_longer(
#				cols = !(sample_id),
#				names_to = "mc",
#				values_to = "mc_score"
#			) |>
#			dplyr::group_by(sample_id) |>
#			dplyr::slice_max(mc_score, n = 6) 
#
			#and then filter by input$selectSample
#	})
}

shinyApp(ui, server) 

