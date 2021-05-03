# Satisfaction survey of the OpenWater 17 in 2020 for the Crowd4SDG project

# Needed Libraries----
library(tidyverse)
library(readxl)
library(gridExtra)
library(gcookbook)
library(labelled)
library(survey)
library(hrbrthemes)
library(viridis)
library (gtsummary)
library(patchwork)
library(hrbrthemes)
library(circlize)
library(networkD3)
# Colors for plots and graphs ----
mycolours1 <- c(
  "Not helpful at all"="#d7191c",
  "Unhelpful"="#fdae61",
  "Neither helpful nor unhelpful"="#ffffbf",
  "Helpful"="#abd9e9",
  "Very helpful"="#2c7bb6",
  "Did not attend & did not watch the playback"="black")
mycolours2 <- c(
  "Strongly disagree"="#d7191c",
  "Disagree"="#fdae61",
  "Neither disagree or agree"="#ffffbf",
  "Agree"="#abd9e9",
  "Strongly agree"="#2c7bb6")

# WORKING ON THE OPEN17 SATISFACTION SURVEY DATA  ----
O17 <- read_excel("d:/Users/masse/Desktop/LAB Marc/Crowd4SDGs/data/Open17Water Challenge Participant Survey.xlsx", sheet = "responses")

# Renaming variables ----
O17r <- O17
names(O17r) <- c("time",
                 "gender",
                 "age",
                 "country_origin",
                 "country_resid",
                 "study_level",
                 "study_field",
                 "start_awareness",
                 "start_decision",
                 "start_process",
                 "team",
                 "team_decision",
                 "team_meet_online",
                 "team_meet_offline",
                 "team_comm_email",
                 "team_comm_whatsapp",
                 "team_comm_slack",
                 "team_comm_goodwall",
                 "team_comm_zoom",
                 "team_comm_wechat",
                 "team_comm_facebook",
                 "team_comm_other_videocall",
                 "team_comm_other_socialmedia",
                 "team_comm_other_app",
                 "team_comm_other_descrip",
                 "team_difficulties",
                 "team_connect_mentors",
                 "team_connect_coordinators",
                 "team_connect_experts",
                 "team_feedback",
                 "team_feedback_satisf",
                 "material_use_handbook",
                 "material_use_curriculum",
                 "material_use_google_folder",
                 "material_use_google_sdginprogress",
                 "issue_zoom",
                 "issue_drive",
                 "issue_goodwall",
                 "issue_slack",
                 "issue_sdginprogress",
                 "issue_internet_open",
                 "rate_method_problemdef",
                 "rate_method_theoryofchange",
                 "rate_method_personas",
                 "rate_method_elevatorpitch",
                 "rate_tool_projectbuilder",
                 "rate_tool_sdginprogress",
                 "rate_tool_kobo",
                 "rate_tool_decidim",
                 "rate_expert_intro",
                 "rate_expert_youngwater",
                 "rate_expert_wmo",
                 "rate_expert_unosat",
                 "rate_expert_dam",
                 "satisf_utility",
                 "satisf_understand_sdg",
                 "satisf_understand_cs",
                 "satisf_understand_crowdsourcing",
                 "satisf_ability_pitch",
                 "satisf_ability_innovate",
                 "satisf_utility_project",
                 "satisf_utility_project",
                 "open_favourite",
                 "open_change",
                 "open_comments")
# Work around demography and profile of participants ----
demo <- subset(O17r, select = c("gender",
                                "age",
                                "country_origin",
                                "country_resid",
                                "study_level",
                                "study_field"))
## Recoding variable
demo$gender<- as.factor(demo$gender)
demo$country_origin <- fct_recode(demo$country_origin,
                                  "Ivory Coast" = "Cote d'Ivoire",
                                  "Italy" = "italy",
                                  "Myanmar" = "Myanmar (Burma)",
                                  "United States" = "The United States"
)
demo$country_resid <- fct_recode(demo$country_resid,
                                 "Italy" = "italy",
                                 "Myanmar" = "Myanmar (Burma)",
                                 "Philippines" = "Quezon City",
                                 "China" = "Shanghai",
                                 "United States" = "United States of America",
                                 "United States" = "USA",
                                 "Myanmar" = "Yangon")
## Recoding demo$study_level
demo$study_level <- fct_recode(demo$study_level,
                               "Undergraduate" = "Am done with Bachelor degree hope I can start master",
                               "Undergraduate" = "Awaiting admission for Masters",
                               "Undergraduate" = "Awaiting Admission for my Masters",
                               "High school" = "Highschool",
                               "Post graduate" = "Post graduate isolated disciplines",
                               "Graduate" = "Recently finished Masters (2020)",
                               "Undergraduate" = "University (Bachelors)",
                               "Graduate" = "University (Masters)",
                               NULL = "Working"
)
demo$study_level <- fct_relevel(
  demo$study_level,
  "High school", "Undergraduate", "Graduate", "Post graduate"
)
demo$study_field <- fct_recode(demo$study_field,
                               "Computer science" = "Computer sciences",
                               "Environmental Science" = "Environment Science",
                               "Environmental science" = "Environmental Geography",
                               "Environmental science" = "Geology/Geochemistry",
                               "Engineering" = "Hydrology and Water Resources Engineering",
                               "Life science" = "Life science (biology, medicine, etc)",
                               "Environmental science" = "Natural resources sciences",
                               "Public administration" = "Public Administration",
                               "Social science" = "Social sciences (law, ethics, sociology, anthropology, economics, etc)",
                               "English studies" = "Teaching English as Foreign Language"
)
var_label(demo) <- list(
  gender = "Gender",
  age = "Age",
  country_origin = "Country of origin",
  country_resid = "Country of residence",
  study_level = "Study level",
  study_field = "Study field")

## Demographics: summary table

demo %>%
  select(gender, age) %>%
  tbl_summary(
  statistic = list(all_continuous() ~ "{median} ({min}, {max})",
                   all_categorical() ~ "{n} ({p}%)"))


demo %>%
  select(study_level) %>%
  tbl_summary(
    statistic = list(all_continuous() ~ "{median} ({min}, {max})",
                     all_categorical() ~ "{n} ({p}%)"))
demo %>%
  select(country_origin) %>%
  tbl_summary(
    statistic = list(all_continuous() ~ "{median} ({min}, {max})",
                     all_categorical() ~ "{n} ({p}%)"))

### I cannot save the plots: 
### as_gt() %>%
    # gt::gtsave(filename = "demo_country_origin.png")

### ggsave("demo_country_origin.png", width = 6, height = 10) 

demo %>%
  select(country_resid) %>%
  tbl_summary(
    statistic = list(all_continuous() ~ "{median} ({min}, {max})",
                     all_categorical() ~ "{n} ({p}%)"))



## Demographics plots ----

### age and gender

demo %>%
  ggplot( aes(x=gender, y=age, fill=gender)) +
  geom_boxplot()+
  scale_fill_brewer()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  scale_y_continuous(breaks = c(10, 15, 20, 25, 30, 35, 40, 45))+
  labs(x = "", y="age", fill = "gender")

### Country of origin

country1 <- demo %>% 
  count(country_origin)

### Compute percentages
country1$fraction = country1$n / sum(country1$n)

### Extract if want to work on excel instead
write.table(country1, file = "d:\\Users\\masse\\Desktop\\Projets\\R-projects\\Crowd4SDG\\country1.csv", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
write.csv(country1, file = "d:\\Users\\masse\\Desktop\\Projets\\R-projects\\Crowd4SDG\\country1.csv")


### Country of origin

country2 <- demo %>% 
  count(country_resid)

### Compute percentages
country2$fraction = country2$n / sum(country2$n)

### Extract if want to work on excel instead
write.table(country2, file = "d:\\Users\\masse\\Desktop\\Projets\\R-projects\\Crowd4SDG\\country2.csv", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
write.csv(country2, file = "d:\\Users\\masse\\Desktop\\Projets\\R-projects\\Crowd4SDG\\country2.csv")



## Plots for country1

### Compute the cumulative percentages (top of each rectangle)
country1$ymax = cumsum(country1$fraction)

### Compute the bottom of each rectangle
country1$ymin = c(0, head(country1$ymax, n=-1))

### Compute label position
country1$labelPosition <- (country1$ymax + country1$ymin) / 2

### Compute a good label
country1$label <- paste0(country1$country_origin,": ", country1$n)

### Make the plot -- and gave up because it is ugly
ggplot(country1, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=country_origin)) +
  geom_rect() +
  geom_text( x=5, aes(y=labelPosition, label=label, color=country_origin), size=2.5) + # x here controls label position (inner / outer)+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")

# Primary work around 017 satisfaction ----

## Let's format the data for 017 satisfaction, starting with the "rate" section where participants had the choice to rank satisfaction from very unsuful to very useuful

O17_satisf <- subset(O17r, select = c("rate_method_problemdef",
                                      "rate_method_theoryofchange",
                                      "rate_method_personas",
                                      "rate_method_elevatorpitch",
                                      "rate_tool_projectbuilder",
                                      "rate_tool_sdginprogress",
                                      "rate_tool_kobo",
                                      "rate_tool_decidim",
                                      "rate_expert_intro",
                                      "rate_expert_youngwater",
                                      "rate_expert_wmo",
                                      "rate_expert_unosat",
                                      "rate_expert_dam",
                                      "satisf_utility",
                                      "satisf_understand_sdg",
                                      "satisf_understand_cs",
                                      "satisf_understand_crowdsourcing",
                                      "satisf_ability_pitch",
                                      "satisf_ability_innovate",
                                      "satisf_utility_project",
                                      "satisf_utility_project",
                                      "open_favourite",
                                      "open_change",
                                      "open_comments"))
                     
rate_cols = grep("rate", colnames(O17_satisf)) #Columns with "Rate" in their title 
O17_rate = O17_satisf[,rate_cols] #Subsetting rate columns
rate = as.data.frame(t(O17_rate)) #Transposes the DF to Matrix
rate$name = rownames(rate) #Reconverts to DF
rate_ok = reshape2::melt(rate, id.vars = "name") #Reshapes by fixing the variable row
rate_ok = rate_ok[,c(1,3)] #Second row is unnecessary
colnames(rate_ok) = c("name", "level") #renaming to avoid ambiguity
rate_ok = reshape2::melt(table(rate_ok)) #Table counts the frequency and then reshapes again
rate_ok$level <- factor(rate_ok$level,levels = c("Not helpful at all", "Unhelpful", "Neither helpful nor unhelpful", "Helpful", "Very helpful","Did not attend & did not watch the playback" ))

### Renaming labels for the legend using facets
rate_method_label <- c("Problem definition", "Theory of change", "Personas", "Elevator Pitch")
names(rate_method_label) <- c("rate_method_problemdef","rate_method_theoryofchange","rate_method_personas","rate_method_elevatorpitch")

rate_tool_label <- c("CS Project Builder", "SDG in progress", "Kobo App", "Decidim")
names(rate_tool_label) <- c("rate_tool_projectbuilder","rate_tool_sdginprogress","rate_tool_kobo","rate_tool_decidim")

rate_expert_label <- c("Urban Water resilience","Young Water Solutions", "Flood Rapid Mapping (UNOSAT)", "Dam Plastics (Case study)")
names(rate_expert_label) <- c("rate_expert_intro","rate_expert_wmo","rate_expert_unosat","rate_expert_dam")
  
# Descriptive analysis: O17 Satisfaction  -----

## Let's have a look at the satisfaction of the workshop's methods

ggplot(subset(rate_ok, name %in% c("rate_method_problemdef",
                                   "rate_method_theoryofchange",
                                   "rate_method_personas",
                                   "rate_method_elevatorpitch")))+ 
  aes(x = level, y = value, fill = level)+
  geom_bar( stat = "identity", position = position_dodge())+
  scale_fill_manual (values = mycolours1)+
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30))+
  scale_x_discrete(labels = NULL, breaks = NULL) +
  labs(x = "", y="Count", fill = "satisfaction level")+
  theme(legend.position = "bottom")+
  geom_text(aes(label= value), position=position_dodge(width=0.9), size = 3, vjust=-0.25)+
  facet_grid(cols = vars(name), labeller = labeller(name = rate_method_label))
  
## Let's have a look at the satisfaction of the workshop's tools

ggplot(subset(rate_ok, name %in% c("rate_tool_projectbuilder","rate_tool_sdginprogress","rate_tool_kobo","rate_tool_decidim")))+
  aes(x = level, y = value, fill = level)+
  geom_bar( stat = "identity", position = position_dodge())+
  scale_fill_manual (values = mycolours1)+
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30))+
  scale_x_discrete(labels = NULL, breaks = NULL) +
  labs(x = "", y="Count", fill = "satisfaction level")+
  theme(legend.position = "bottom")+
  geom_text(aes(label= value), position=position_dodge(width=0.9), size = 3, vjust=-0.25)+
  facet_grid(cols = vars(name), labeller = labeller(name = rate_tool_label))

## Let's have a look at the satisfaction of the workshop's presentation of experts

ggplot(subset(rate_ok, name %in% c("rate_expert_intro","rate_expert_wmo","rate_expert_unosat","rate_expert_dam")))+
  aes(x = level, y = value, fill = level)+
  geom_bar( stat = "identity", position = position_dodge())+
  scale_fill_manual (values = mycolours1)+
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30))+
  scale_x_discrete(labels = NULL, breaks = NULL) +
  labs(x = "", y="Count", fill = "satisfaction level")+
  theme(legend.position = "bottom", strip.text = element_text(size=9,lineheight=10))+
  geom_text(aes(label= value), position=position_dodge(width=0.9), size = 3, vjust=-0.25)+
  facet_grid(cols = vars(name), labeller = labeller(name = rate_expert_label))


## Formating the data for 017 satisfaction focusing on the "satisf_" section where participants had the choice to rank satisfaction  and learning using Likert scale ranging from strongly disagree to strongly agree

satisf_cols = grep("satisf", colnames(O17_satisf)) #Columns with "Satisf" in their title 
O17_satisf = O17_satisf[,satisf_cols] #Subsetting satisf columns
satisf = as.data.frame(t(O17_satisf)) #Transposes the DF to Matrix
satisf$name_satisf = rownames(satisf) #Reconverts to DF
satisf_ok = reshape2::melt(satisf, id.vars = "name_satisf") #Reshapes by fixing the variable row
satisf_ok = satisf_ok[,c(1,3)] #Second row is unnecessary
colnames(satisf_ok) = c("name_satisf", "level_satisf") #renaming to avoid ambiguity
satisf_ok = reshape2::melt(table(satisf_ok)) #Table counts the frequency and then reshapes again



satisf_ok$level_satisf <- factor(satisf_ok$level_satisf,levels = c("Strongly disagree",
                                                 "Disagree",
                                                 "Neither disagree or agree",
                                                 "Agree",
                                                 "Strongly agree"))

## ## Let's have a look at the global perceived satisfaction "O17 met the expectations"

ggplot(subset(satisf_ok, name_satisf %in% c("satisf_utility")))+
  aes(x = level_satisf, y = value, fill = level_satisf)+
  geom_bar( stat = "identity", position = position_dodge())+
  scale_fill_manual (values = mycolours2)+
  scale_y_continuous(limits = c(0, 28), breaks = c(0, 5, 10, 15, 20, 25, 30))+
  scale_x_discrete(labels = NULL, breaks = NULL) +
  labs(x = "'The #Open17Water Challenge met the expectations I had of the program.'",fill = " ")+
  theme(legend.position = "bottom")+
  geom_text(aes(label= value), position=position_dodge(width=0.9), size = 3, vjust=-0.25)

## Let's have a look at the perceived learning - aka increased in understanding


satisf_understanding_label <- c("the SDGs","Citizen Science", "Crowdsourcing") # Creating labels for the legend of facets
names(satisf_understanding_label) <- c("satisf_understand_sdg","satisf_understand_cs","satisf_understand_crowdsourcing")



ggplot(subset(satisf_ok, name_satisf %in% c("satisf_understand_sdg",
                                            "satisf_understand_cs",
                                            "satisf_understand_crowdsourcing")))+
  aes(x = level_satisf, y = value, fill = level_satisf)+
  geom_bar( stat = "identity", position = position_dodge())+
  scale_fill_manual (values = mycolours2)+
  scale_y_continuous(limits = c(0, 28), breaks = c(0, 5, 10, 15, 20, 25, 30))+
  scale_x_discrete(labels = NULL, breaks = NULL) +
  labs(x = "'Since participating in the #Open17Water Challenge, my understanding of ______ has increased.'",fill = " ")+
  theme(legend.position = "bottom")+
  geom_text(aes(label= value), position=position_dodge(width=0.9), size = 3, vjust=-0.25)+
  facet_grid(cols = vars(name_satisf), labeller = labeller(name_satisf = satisf_understanding_label))


## Let's have a look at the perceived capacity building - aka increased ability


satisf_ability_label <- c("Pitch","Innovate") # Creating labels for the legend of facets
names(satisf_ability_label) <- c("satisf_ability_pitch", "satisf_ability_innovate")


ggplot(subset(satisf_ok, name_satisf %in% c("satisf_ability_pitch", "satisf_ability_innovate")))+
  aes(x = level_satisf, y = value, fill = level_satisf)+
  geom_bar( stat = "identity", position = position_dodge())+
  scale_fill_manual (values = mycolours2)+
  scale_y_continuous(limits = c(0, 28), breaks = c(0, 5, 10, 15, 20, 25, 30))+
  scale_x_discrete(labels = NULL, breaks = NULL) +
  labs(x = "'Since participating in the #Open17Water Challenge, my ability to ______ has increased.'",fill = " ")+
  theme(legend.position = "bottom")+
  geom_text(aes(label= value), position=position_dodge(width=0.9), size = 3, vjust=-0.25)+
  facet_grid(cols = vars(name_satisf), labeller = labeller(name_satisf = satisf_ability_label))

# Now let's look at how the perceived gained knowledge will allow participants to develop their projects

ggplot(subset(satisf_ok, name_satisf %in% c("satisf_utility_project")))+
  aes(x = level_satisf, y = value, fill = level_satisf)+
  geom_bar( stat = "identity", position = position_dodge())+
  scale_fill_manual (values = mycolours2)+
  scale_y_continuous(limits = c(0, 28), breaks = c(0, 5, 10, 15, 20, 25, 30))+
  scale_x_discrete(labels = NULL, breaks = NULL) +
  labs(x = "'I'm going to use the knowledge that I have gained during the #Open17Water Challenge to further develop my project.'",fill = " ")+
  theme(legend.position = "bottom")+
  geom_text(aes(label= value), position=position_dodge(width=0.9), size = 3, vjust=-0.25)

## Exporting data and plots

write.table(demo, file = "d:\\Users\\masse\\Desktop\\Projets\\R-projects\\Crowd4SDG\\O17_demographics.csv", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
write.csv(demo, file = "d:\\Users\\masse\\Desktop\\Projets\\R-projects\\Crowd4SDG\\O17_demographics.csv")


# WORKING ON OPEN17 DATA FROM AMUDHA - PARTICIPANTS AND TEAMS -----
#### This is different from the survey data

participants <- read_csv("Open17_participants.csv")


## Recoding variables

participants$team <- as.factor(participants$team)
participants$gender <- as.factor(participants$gender)
participants$country_origine <- as.factor(participants$country_origine)
participants$country_resid <- as.factor(participants$country_resid)


participants <- subset(participants, select = c("team",
                                        "gender",
                                        "age",
                                        "country_origine",
                                        "country_resid"))

participants$country_resid <- fct_recode(participants$country_resid,
                                "China" = "Shanghai")


## Univariate analysis

var_label(participants) <- list(
  gender = "Gender",
  age = "Age",
  country_origine = "Country of origin",
  country_resid = "Country of residence")

count_country_o <- participants %>% count(country_origine)
names(count_country_o) <- c("co", "countco")

count_country_r <- participants %>% count(country_resid)
names(count_country_r) <- c("cr", "countcr")

count_team <- participants %>% count(team)
names(count_team)<- c("team", "countteam")



## Demographics: summary table

participants %>%
  select(gender, age) %>%
  tbl_summary(
    statistic = list(all_continuous() ~ "{median} ({min}, {max})",
                     all_categorical() ~ "{n} ({p}%)"))

ggplot(count_country_o)+
  aes(x = reorder(co,countco), y= countco)+
  geom_bar(stat="identity",width = 0.9, position = position_dodge())+
  scale_y_continuous(limits = c(0, 15), breaks = c(0, 5, 10, 15))+
  labs(y = "Count",x = "Country of origin")+
  theme_bw() +
  geom_text(aes(label = countco),
            position = position_dodge(width = 0.9), size = 2.5, hjust = -0.5)+
  coord_flip()
  

ggplot(count_country_r)+
  aes(x = reorder(cr,countcr), y= countcr)+
  geom_bar(stat="identity",width = 0.9, position = position_dodge())+
  scale_y_continuous(limits = c(0, 15), breaks = c(0, 5, 10, 15))+
  labs(y = "Count",x = "Country of residence")+
  theme_bw() +
  geom_text(aes(label = countcr),
            position = position_dodge(width = 0.9), size = 2.5, hjust = -0.5)+
  coord_flip()


ggplot(count_team)+
  aes(x = reorder(team,countteam), y= countteam)+
  geom_bar(stat="identity",width = 0.9, position = position_dodge())+
  scale_y_continuous(limits = c(0, 5), breaks = c(0,1,2,3,4,5))+
  labs(y = "Count",x = "", title = "Team size")+
  theme_bw()+
  coord_flip()

# p <- wilcox.test(x,y)$p.value
# title(paste('p =', signif(p,2))

participants %>%
  ggplot( aes(x=gender, y=age, fill=gender)) +
  geom_boxplot()+
  scale_fill_brewer()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  scale_y_continuous(breaks = c(10, 15, 20, 25, 30, 35, 40, 45))+
  labs(x = "", y="age", fill = "gender")+
  theme_bw()


