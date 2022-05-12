# Experimental Script for Dashboard Update
# Packages ----
library(shiny)  # Required to run any Shiny app
library(tidyverse) # read_csv
library(ggplot2)  # For creating pretty plots
library(lubridate) # dates
library(plyr)
library(zoo)
library(reshape)
library(plotly)

#Loading-data
data_pre <- read.csv("https://raw.githubusercontent.com/akvariants/modeldata/main/metadata.alaska.recent.csv", header=T) #variant data pango version 3.1.17
variant_def <- read.csv("https://raw.githubusercontent.com/akvariants/modeldata/main/variant_def.csv", header=T) #variant classification
data <- merge(data_pre,variant_def, by.x = "lineage", by.y = "Pangolin", all.x = TRUE) 
data$date <- as.Date(data$`Collection.date`, "%m/%d/%y")

time_summary <- subset(data, data$date >= as.Date("2022-01-01"))
time_summary$lineage_cor <- time_summary$lineage 
time_summary$lineage_cor[time_summary$lineage_cor=="BA.2.1"] <- "BA.2"
time_summary$lineage_cor[time_summary$lineage_cor=="BA.2.3"] <- "BA.2"
time_summary$lineage_cor[time_summary$lineage_cor=="BA.2.3.1"] <- "BA.2"
time_summary$lineage_cor[time_summary$lineage_cor=="BA.2.3.2"] <- "BA.2"
time_summary$lineage_cor[time_summary$lineage_cor=="BA.2.7"] <- "BA.2"
time_summary$lineage_cor[time_summary$lineage_cor=="BA.2.9"] <- "BA.2"
time_summary$lineage_cor[time_summary$lineage_cor=="BA.2.10"] <- "BA.2"
time_summary$lineage_cor[time_summary$lineage_cor=="BA.2.12"] <- "BA.2"
time_summary$lineage_cor[time_summary$lineage_cor=="BA.2.12.1"] <- "BA.2"


time_plot_all <- as.data.frame.matrix(table(time_summary$date,time_summary$lineage_cor))
time_plot_relative_all <- as.data.frame.matrix(t(apply(time_plot_all, 1, function(x) x/sum(x))))
time_plot_all$date <- rownames(time_plot_all)
time_plot_all <- gather(time_plot_all, key, value, -date)
time_plot_relative_all$date <- rownames(time_plot_relative_all)
time_plot_relative_all <- gather(time_plot_relative_all, key, value, -date)
time_plot_relative_all$count <- time_plot_all$value

Omicron_cases_freq <- subset(time_plot_relative_all, time_plot_relative_all$key =="BA.2" &
                               time_plot_relative_all$date >= as.Date("2022-01-09") &
                               time_plot_relative_all$date <= as.Date("2022-04-16"))

Omicron_cases_freq$date <- as.Date(Omicron_cases_freq$date)
omicron_logit <- glm(value~date, family=binomial(link ='logit'), data=Omicron_cases_freq)
Omicron_cases_freq <- cbind(Omicron_cases_freq,as.data.frame(omicron_logit$fitted.values))
colnames(Omicron_cases_freq) <- c("date","key","Proportion of sequenced cases","count","fitted")
Omicron_cases_freq <- add_column(Omicron_cases_freq, mod_se = predict(omicron_logit, newdata = Omicron_cases_freq, type = 'response',
                                                                      se.fit = TRUE)$se.fit)
Omicron_cases_freq_most_recent <- tail(Omicron_cases_freq, n=1)
model_plot <- ggplot(Omicron_cases_freq, aes(x=date, y =`Proportion of sequenced cases`)) +
  geom_point(alpha=0.2, aes(color="Daily BA.2% of Sequenced Cases"))+
  geom_smooth(method = "glm", method.args = list(family=binomial(link = 'logit')),formula=y~x, fill="#b07aa1",aes(color="Predicted Prevalence"))+
  xlab("Date")+
  ylab("Percent of Sequenced Cases")+
  theme_minimal()+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  scale_color_manual(name = "",values =c("#b07aa1","#b07aa1")) 
#    ggtitle("Estimated prevalence of BA.2 from first detection") +
#    theme(plot.title = element_text(color="#00798c",size=14, face="bold",hjust = 0.5))


fitted <- as.data.frame(fitted(omicron_logit))
AKcases <- read.csv("https://raw.githubusercontent.com/akvariants/modeldata/main/AKcases.recent.csv", header= T)
AKcases$Onset_Date <- as.Date(AKcases$Onset_Date, "%m/%d/%y")
cases <- as.data.frame(count(AKcases$Onset_Date))
colnames(cases) <- c("onset_date","new_confirmed_cases")
cases$onset_date <- as.Date(cases$onset_date, "%m/%d/%y")
cases$`Reported Cases` <- rollmean(cases$new_confirmed_cases, k = 7, fill = NA)
Ba2_cases <- subset(cases,cases$onset_date >= as.Date("2022-01-09") &
                      cases$onset_date <= as.Date("2022-04-16"))
Ba2_cases$`Cases of BA.2` <- fitted$`fitted(omicron_logit)`* Ba2_cases$`Reported Cases`
Ba2_cases_melt <- melt(Ba2_cases[,c(1,3,4)],id="onset_date")
colnames(Ba2_cases_melt) <- c("Onset Date", "Type", "Cases")

cases_plot <- ggplot(Ba2_cases_melt, aes(x = `Onset Date`, y = Cases, fill=Type)) +
  geom_area(position="identity")+
  ylab("Cases") +
  xlab("Onset Date")+
  theme_minimal() +
  scale_fill_manual(name = "COVID-19 cases", labels = c("Reported Cases", "Cases of BA.2"),values =c("#D3D3D3", "#de97ca"))
#ggtitle("Estimated Cases Attributed to BA.2")+
#    theme(plot.margin =  unit(c(1,1,0,1), "lines"),legend.position = c(0.8, 0.9),
#          plot.title = element_text(color="#00798c",size=16, face="bold",hjust = 0.5),
#          text =element_text(size=16))


# ui.R ----
ui <- fluidPage(
  titlePanel(h2("Estimating Prevalence of Variants", align='center')),
  headerPanel(h3("Here you can find information on the current emerging lineage in Alaska and its prevalence estimated using logistic modelling.", align='center')),
  tags$head(tags$style('h4 {color:#00798c;}')),
  mainPanel(fluidRow(splitLayout(cellWidths = c("70%", "30%"),
                                 h4("BA.2 Estimated Prevalence from 1st Detection"), h4("Current Prevalence of BA.2"))),
            fluidRow(splitLayout(cellWidths = c("70%", "30%"),
                                 plotlyOutput("plot"), plotOutput("plot3"))),
            h4("Estimated Cases Attributed to BA.2", align='center'),
            plotlyOutput("plot2"), align='center', width=12)
  
)

# server.R ----
server <- function(input, output) {
  output$plot <- renderPlotly(
    ggplotly(model_plot) %>%
      layout(legend = list(orientation = "v", x = 0.2, y = 1))
  )
  output$plot2 <- renderPlotly(
    ggplotly(cases_plot)%>%
      layout(legend = list(orientation = "v", x = 0.5, y = 1))
  )
  output$plot3 <- renderPlot(
    ggplot(Omicron_cases_freq_most_recent,aes(x=as.Date(date)))+
      geom_boxplot(aes(lower=fitted-1.96*mod_se,upper=fitted+1.96*mod_se,middle=fitted,ymin=fitted-3*mod_se,ymax=fitted+3*mod_se),stat="identity",
                   fill="#e1cddc", color="#b07aa1") +
      ylab("Estimated Prevalence of BA.2") +
      xlab("Date")+
      theme_minimal()+
      scale_y_continuous(labels = scales::percent_format(accuracy = 1),limits=c(0,1.05)) +
      scale_x_date(date_labels = "%b %d %y",expand = c(0,0))+
      theme(plot.margin =  unit(c(1,1,0,1), "lines"),legend.position = c(0.8, 0.9),
            text =element_text(size=16))
  )
}


# Run the app ----
shinyApp(ui = ui, server = server)