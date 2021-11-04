library(shiny)
library(ggplot2)
library(tidyverse)
####################

ui <- fluidPage(
  #Create input for domain (single) and variable (multiple)
  selectInput("measure", "Measure", choices = ""),
  selectInput("varSelection", "Variables", multiple = T, choices = ""),
  
  #Set plot for output
  plotOutput("myPlot")
)

server <- function(input, output, session) {
  S1 <- read_csv("S10.42021-10-30.csv")
  S2 <- read_csv("S20.42021-10-30.csv")
  S3 <- read_csv("S30.42021-10-30.csv")
  S4 <- read_csv("S40.42021-10-30.csv")
  SS=rbind(S1,S2,S3,S4)
  
  SS$p=as.factor(SS$p)
  
  SS=SS%>%mutate(Measure=X1)
  SS=SS%>%select(-X1)
  SS=SS%>%select(-n,-pr)
  #Load your data
  #Update the domain input according to the data
  updateSelectInput(session, "measure", choices = unique(SS$Measure) )
  
  #Update the variable list (assumed all but d and year are variables of interest)
  updateSelectInput(session, "varSelection", 
                    choices = colnames(SS %>% select(-Case, -N,-p,-c,-Measure)))
  
  #Load the chart function
  draw_chart <- function(SS, listv, d){

    SS2=SS%>%filter(Measure==d)
    
    SS2=SS2%>%select(-Measure,-c,-ATT.true)
    
    SS2=SS2%>%gather(key="Method",value = d,-Case,-N,-p)
   
     R=range(SS2$d)
    
    SS2 <- SS2 %>%
      filter(Method %in% listv)
  
    
    SS2%>%ggplot(aes(x=N,y=d,color=Method))+
      geom_line()+ylim(R)+
      facet_grid(p~Case)+ggtitle("")

  }
 
   #Render the plot
  output$myPlot = renderPlot({
    draw_chart(SS, input$varSelection, input$measure)
  })

  }


shinyApp(ui, server)

