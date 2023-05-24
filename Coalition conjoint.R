#####Prepare workspace#####
rm(list=ls())
library(tidyverse)
library(rio)
library(cregg)
library(marginaleffects)
library(scales)
library(brms)

theme_plots <- function() {
  theme_minimal(base_family = "IBM Plex Sans Condensed") +
    theme(panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"),
          strip.background = element_rect(fill = "grey95", color = NA),
          legend.title = element_text(face = "bold"))
}

update_geom_defaults("label", 
                     list(family = "IBM Plex Sans Condensed"))
update_geom_defaults(ggtext::GeomRichText, 
                     list(family = "IBM Plex Sans Condensed"))


#####Import and clean data#####
dta <- import('Badanie CAWI-conjoint_kwiecień2023-zbiór.sav')

dta <- dta %>%
  mutate(vote = as.factor(case_when(
    `Q1` == 1 ~ "Yes",
    `Q1` == 2 ~ "Yes",
    `Q1` == 3 ~ "No",
    `Q1` == 4 ~ "No")),
    votefor = as.factor(case_when(
      `Q2` == 1 ~ "KO",
      `Q2` == 2 ~ "PiS",
      `Q2` == 3 ~ "Lewica",
      `Q2` == 4 ~ "Konfederacja",
      `Q2` == 5 ~ "Polska 2050",
      `Q2` == 6 ~ "PSL-KP",
      `Q2` == 88 ~ "Other",
      `Q2` == 98 ~ "Don't know")),
    ubi = coalesce(Q4A, Q4B),
    ubi = replace(ubi, ubi %in% 98:99, NA),
    ubi = as.ordered(ubi),
    ubitype = as.factor(case_when(
      `WERSJA_LOS` == 1 ~ "Neutral wording",
      `WERSJA_LOS` == 2 ~ "Nativist wording")),
    goodbad = replace(Q10, Q10 %in% 98:99, NA),
    goodbad = as.ordered(goodbad)
  )


#####UBI analysis#####
bprior <- c(prior_string("normal(0,1)", class="b"))

ubi_mod <- brm(ubi | weights(waga) ~ votefor*ubitype, data=dta, family=cumulative("logit"), prior=bprior, sample_prior=TRUE, backend="cmdstanr", chains=3, cores=3, threads = threading(3))

my_labels <- c('Strongly\nsupport', 'Support', 'Neither support\nnor oppose', 'Oppose', 'Strongly oppose')

ubi_pred <- predictions(ubi_mod, newdata=datagrid(votefor=c("PiS", "KO", "Konfederacja", "Lewica", "Polska 2050", "PSL-KP"), ubitype=c("Neutral wording", "Nativist wording")), conf_level = 0.8) %>%
  #posterior_draws() %>%
  ggplot() +
  geom_pointrange(aes(x=group, y=estimate, ymin=conf.low, ymax=conf.high, color=ubitype), position=position_dodge(width=0.5)) +
  scale_color_manual(values=c("red", "blue")) +
  scale_x_discrete(labels=my_labels) +
  facet_wrap(~votefor) +
  labs(x="", y="Predicted probability", color="") +
  theme_plots()
ggsave(ubi_pred, file = "ubi_pred.png", width = 8, height = 5, units = "cm", dpi = 320, scale = 4.5, bg="white")

ubi_pred_overall <- predictions(ubi_mod, newdata=datagrid(votefor=c("PiS", "KO", "Konfederacja", "Lewica", "Polska 2050", "PSL-KP")), conf_level = 0.8) %>%
  #posterior_draws() %>%
  ggplot() +
  geom_pointrange(aes(x=group, y=estimate, ymin=conf.low, ymax=conf.high, color=ubitype), position=position_dodge(width=0.5), color="purple") +
  scale_x_discrete(labels=my_labels) +
  facet_wrap(~votefor) +
  labs(x="", y="Predicted probability", color="") +
  theme_plots()
ggsave(ubi_pred_overall, file = "ubi_pred_overall.png", width = 8, height = 5, units = "cm", dpi = 320, scale = 4.5, bg="white")

ubi_comp <- comparisons(ubi_mod, variables=list("ubitype"=c("Neutral wording", "Nativist wording")), by="votefor", newdata=datagrid(votefor=c("PiS", "KO", "Konfederacja", "Lewica", "Polska 2050", "PSL-KP")), conf_level = 0.8) %>%
  ggplot() +
  geom_hline(aes(yintercept=0), colour="gray40", linetype="dotted") +
  geom_pointrange(aes(x=group, y=estimate, ymin=conf.low, ymax=conf.high), position=position_dodge(width=0.5), color="blue") +
  scale_x_discrete(labels=my_labels) +
  facet_wrap(~votefor) +
  labs(x="", y="Contrast", color="") +
  theme_plots()
ggsave(ubi_comp, file = "ubi_comp.png", width = 8, height = 5, units = "cm", dpi = 320, scale = 4.5, bg="white")

ubi_support <- ggplot(data=subset(dta, !is.na(ubi)), aes(x=ubi, fill=ubi)) + 
  geom_bar(aes(y = (..count..)/sum(..count..))) + 
  scale_y_continuous(labels = percent) +
  scale_x_discrete(labels=my_labels) +
  scale_fill_brewer(palette="RdYlGn", guide=FALSE, direction=-1) + 
  labs(x="Attitude to UBI", y="Level of support", color="") +
  theme_plots()
ggsave(ubi_support, file = "ubi_support.png", width = 8, height = 5, units = "cm", dpi = 320, scale = 4.5, bg="white")


#####Conjoint analysis#####
list1 <- list(
  vote = list(~ CJ_1_1_a + CJ_2_1_a + CJ_3_1_a + CJ_4_1_a + CJ_5_1_a + CJ_6_1_a + CJ_7_1_a + CJ_8_1_a + CJ_9_1_a + CJ_10_1_a + + CJ_11_1_a + CJ_12_1_a,
                ~ CJ_1_2_a + CJ_2_2_a + CJ_3_2_a + CJ_4_2_a + CJ_5_2_a + CJ_6_2_a + CJ_7_2_a + CJ_8_2_a + CJ_9_2_a + CJ_10_2_a + + CJ_11_2_a + CJ_12_2_a
  ),
  gender = list(~ CJ_1_1_b + CJ_2_1_b + CJ_3_1_b + CJ_4_1_b + CJ_5_1_b + CJ_6_1_b + CJ_7_1_b + CJ_8_1_b + CJ_9_1_b + CJ_10_1_b + + CJ_11_1_b + CJ_12_1_b,
              ~ CJ_1_2_b + CJ_2_2_b + CJ_3_2_b + CJ_4_2_b + CJ_5_2_b + CJ_6_2_b + CJ_7_2_b + CJ_8_2_b + CJ_9_2_b + CJ_10_2_b + + CJ_11_2_b + CJ_12_2_b
  ),
  age = list(~ CJ_1_1_c + CJ_2_1_c + CJ_3_1_c + CJ_4_1_c + CJ_5_1_c + CJ_6_1_c + CJ_7_1_c + CJ_8_1_c + CJ_9_1_c + CJ_10_1_c + + CJ_11_1_c + CJ_12_1_c,
              ~ CJ_1_2_c + CJ_2_2_c + CJ_3_2_c + CJ_4_2_c + CJ_5_2_c + CJ_6_2_c + CJ_7_2_c + CJ_8_2_c + CJ_9_2_c + CJ_10_2_c + + CJ_11_2_c + CJ_12_2_c
  ),
  exper = list(~ CJ_1_1_d + CJ_2_1_d + CJ_3_1_d + CJ_4_1_d + CJ_5_1_d + CJ_6_1_d + CJ_7_1_d + CJ_8_1_d + CJ_9_1_d + CJ_10_1_d + + CJ_11_1_d + CJ_12_1_d,
              ~ CJ_1_2_d + CJ_2_2_d + CJ_3_2_d + CJ_4_2_d + CJ_5_2_d + CJ_6_2_d + CJ_7_2_d + CJ_8_2_d + CJ_9_2_d + CJ_10_2_d + + CJ_11_2_d + CJ_12_2_d
  ),
  listpl = list(~ CJ_1_1_e + CJ_2_1_e + CJ_3_1_e + CJ_4_1_e + CJ_5_1_e + CJ_6_1_e + CJ_7_1_e + CJ_8_1_e + CJ_9_1_e + CJ_10_1_e + + CJ_11_1_e + CJ_12_1_e,
              ~ CJ_1_2_e + CJ_2_2_e + CJ_3_2_e + CJ_4_2_e + CJ_5_2_e + CJ_6_2_e + CJ_7_2_e + CJ_8_2_e + CJ_9_2_e + CJ_10_2_e + + CJ_11_2_e + CJ_12_2_e
  )
)

list2 <- list(choice = ~ Q3_1 + Q3_2 + Q3_3 + Q3_4 + Q3_5 + Q3_6 + Q3_7 + Q3_8 + Q3_9 + Q3_10 + Q3_11 + Q3_12)


dta_long <- cj_tidy(dta,
                    profile_variables = list1,
                    task_variables = list2,
                    id = ~ id)

dta_long$chosen <- ifelse((dta_long$profile == "A" & dta_long$choice == 1) | 
                            (dta_long$profile == "B" & dta_long$choice == 2), 1, 0)

dta_long <- tibble(dta_long) %>%
  mutate(vote = factor(vote, labels=c("PiS", "Konfederacja", "KO", "PSL", "Lewica", "Polska 2050", 
                                        "Koalicja Lewica i Polska 2050", "Koalicja Polska 2050 i PSL",
                                        "Koalicja KO i Lewica", "Koalicja KO i PSL", "Koalicja KO, PSL i Polska 2050",
                                        "Koalicja KO, PSL, Polska 2050 i Lewica")),
         listpl = factor(listpl, labels=c("1", "4", "9")),
         gender = factor(gender, labels=c("Male", "Female")),
         age = factor(age, labels=c("23", "31", "39", "46", "57")),
         exper = factor(exper, labels=c("New candidate", "Current MP"))
  )

mod_amce <- cj(dta_long, chosen ~
            vote + gender + age + exper + listpl, by = ~votefor,
          id = ~pair, weights= ~waga, estimate="amce",
          feature_order = c("vote", "gender", "age", "exper", "listpl"))

plot(mod_amce, size=2) +
  facet_wrap(~votefor, ncol = 3L) +
  theme_plots() 

dta_long_koalicja <- tibble(dta_long) %>%
  mutate(vote = relevel(vote, ref = "Koalicja KO, PSL, Polska 2050 i Lewica")
  )

mod_koalicja_amce <- cj(dta_long_koalicja, chosen ~
                 vote + gender + age + exper + listpl, by = ~votefor,
               id = ~pair, weights= ~waga, estimate="amce",
               feature_order = c("vote", "gender", "age", "exper", "listpl"))

plot(mod_koalicja_amce, size=2) +
  facet_wrap(~votefor, ncol = 3L) +
  theme_plots() 



goodbad_mod <- brm(goodbad | weights(waga) ~ votefor, data=dta, family=cumulative("logit"), prior=bprior, sample_prior=TRUE, backend="cmdstanr", chains=3, cores=3, threads = threading(3))

my_labels <- c('Definitely a\nbad person', 'Probably a\nbad person', 'Cannot say whether\n good or bad', 'Probably a\ngood person', 'Definitely a\ngood person')

goodbad_pred <- predictions(goodbad_mod, newdata=datagrid(votefor=c("PiS", "KO", "Konfederacja", "Lewica", "Polska 2050", "PSL-KP")), conf_level = 0.8) %>%
  ggplot() +
  geom_pointrange(aes(x=group, y=estimate, ymin=conf.low, ymax=conf.high), color="red") +
  scale_x_discrete(labels=my_labels) +
  facet_wrap(~votefor) +
  labs(x="", y="Predicted probability", color="", title="If you hear that someone is a PiS supporter and you have no other information about them, which of the following views would you take of them?")+
  theme_plots()
ggsave(goodbad_pred, file = "goodbad_pred.png", width = 8, height = 5, units = "cm", dpi = 320, scale = 5, bg="white")
