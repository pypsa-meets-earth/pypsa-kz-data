# for the first run use: install.packages(c("tidyverse", "lubridate"))
# to install all the needed libraries

rm(list = ls())

library(tidyverse)
library(lubridate)

data_dir <- "./"

# taken from "pypsa-earth/data/ssp2-2.6/2030/era5_20*"
data_fl_list <- list("Asia_demand_2011.csv", "Asia_demand_2013.csv", 
    "Asia_demand_2018.csv")

# load data --------------------------------------------------------------------
demand_gegis_list <- vector(mode = "list", length = length(data_fl_list))

for ( i in seq(along.with = demand_gegis_list) ){
    data_fl <- data_fl_list[[i]]
    demand_gegis_list[[i]] <- read.table(file.path(data_dir, data_fl), 
        sep = ";", header = TRUE, skipNul = TRUE)
}

demand_gegis_df <- bind_rows(demand_gegis_list)
print(head(demand_gegis_df))


# KZ oficial statistics
demand_korem_2011 <-  data.frame(
    year = 2011,
    month = seq(from = 1, to = 12, by = 1),
    demand_korem = c(8633.4, 7665.5, 8048.1, 6796.6,
        6573.8, 6340.5, 6536.7, 6578.7, 6501.8,
        NA, NA, NA)
) %>%
mutate(
    demand_korem_nz = (demand_korem - mean(demand_korem, na.rm = TRUE))/sd(demand_korem, na.rm = TRUE)
)

demand_korem_2013 <-  data.frame(
    year = 2013,
    month = seq(from = 1, to = 12, by = 1),
    demand_korem = c(8774.3, 7864.5, 8005.9, 7041.5, 
        6894.6, 6520.2, 6798.9, 6820.5, 6733.5,
        7647.1, 7889.0, 8667.2)
) %>%
mutate(
    demand_korem_nz = (demand_korem - mean(demand_korem))/sd(demand_korem)
)
demand_korem_2018 <-  data.frame(
    year = 2018,
    month = seq(from = 1, to = 12, by = 1),
    demand_korem = c(9867.6, 8833.3, 9062.4,
        8105.1, 7956.0, 7606.8,
        8150.9, 7882.1, 7743.3,
        8663.0, 9262.9, 9996.4)
) %>%
mutate(
    demand_korem_nz = (demand_korem - mean(demand_korem))/sd(demand_korem)
)

# cross-check against aggregated values ----------------------------------------

# 2011 ------------------------------------
# 24347,0 for jan-march ("Корем_отчет за 9 мес.doc")
sum(demand_korem_2011[1:3, "demand_korem"])
# 19710,9 for apr-jun ("Корем_отчет за 9 мес.doc")
sum(demand_korem_2011[4:6, "demand_korem"])
# 19617,2 for jul-sep ("Корем_отчет за 9 мес.doc")
sum(demand_korem_2011[7:9, "demand_korem"])
# 63675,1 for jan-sep ("Корем_отчет за 9 мес.doc")
sum(demand_korem_2011[1:9, "demand_korem"])

# 2013 ------------------------------------
# 89640,8 ("Годовой_отчет_2013_АО_КОРЭМ.pdf.pdf")
sum(demand_korem_2013[, "demand_korem"])

# 2018 ------------------------------------
# 103228.3 according to the Annual report ("Годовой_отчет_2018_г.pdf")
sum(demand_korem_2018[, "demand_korem"])


# compose the integral dataset -------------------------------------------------
demand_gegis_kz <- demand_gegis_df

demand_gegis_kz_mth_aggr <- demand_gegis_kz %>%
    dplyr::filter(region_name == "Kazakhstan") %>%     
    mutate(
        year = year(ymd_hms(time)),
        month = month(ymd_hms(time)),
        day = day(ymd_hms(time))
    ) %>%
    group_by(year, month) %>%
    summarize(
        demand_gegis = sum(Electricity.demand)/1e3
    ) %>%
    mutate(
        demand_gegis_nz = (demand_gegis - mean(demand_gegis))/sd(demand_gegis)
    )
print(head(demand_gegis_kz_mth_aggr))

dem_df <- full_join(
    demand_gegis_kz_mth_aggr, 
    bind_rows(
        list(demand_korem_2011, demand_korem_2013, demand_korem_2018)) 
    ) %>%
    mutate(scale = demand_gegis/demand_korem)
print(head(dem_df))

res_df <- dem_df %>%
    select(year, month, demand_gegis, demand_korem)
write.csv(res_df, "kz_demand_validation.csv")


# graphics ---------------------------------------------------------------------
pl <- dem_df %>% 
    ggplot(aes(x = as.integer(month), y = demand_gegis_nz)) +
    facet_grid(rows = vars(year)) +
    geom_line(col = "red") +
    geom_line(aes(y = demand_korem_nz), col = "forestgreen") +
    scale_x_continuous(
        n.breaks = 12, labels = c("", month.abb)
    ) +
    scale_colour_manual(
        breaks = c("a", "b"), 
        values = c("red", "forestgreen"),
        labels = c("Measurements", "Model")
    ) +
    # title("KZ") +
    xlab("") +
    ylab("Normalized electricity demand") +
    theme_bw() +
    theme(text=element_text(size=16)) #+
    # theme(panel.grid.major = element_line(colour = "gray90", size = 0.5))
ggsave("demand_profile_kz.png", width = 7, heigh = 4.5, dpi = 300)

   


