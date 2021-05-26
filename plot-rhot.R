pdf <- data.frame(date = wind$target_end_date, p = x$ρ)
p <- ggplot(pdf, aes(x = date, y = p)) + geom_point() + xlab("Date") + ylab(expression(paste(rho[t], " = Pr (Removal is reported)")))
ggsave("rhot-estimate.pdf", p, width = 7, height = 3)
