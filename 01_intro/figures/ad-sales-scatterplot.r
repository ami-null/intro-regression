library(ggplot2)
ad_data <- data.frame(
      spend = c(1,2,3,4,5,6,7,8,9,10),
      sales = c(1.2,2.3,2.0,3.8,4.0,5.3,5.8,6.9,7.5,8.0)
    )
ggplot(ad_data, aes(x = spend, y = sales)) +
      geom_point(size = 2, color = "#1b9e77") +
      geom_smooth(method = "lm", se = FALSE, color = "#d95f02") +
      labs(x = "Advertising Spend ($1000s)", y = "Sales ($1000s)") +
      theme_minimal() +
      theme(
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background  = element_rect(fill = "transparent", colour = NA)
          )
ggsave("ad_sales_scatter.pdf", bg = "transparent")
