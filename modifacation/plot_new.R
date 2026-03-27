library(ggplot2)
library(cowplot)
library(gridGraphics)
library(grid)
library(ggplotify)
library(tidyverse)
library(readr)

##############
## Figure 1
# (a)
IGRA_1991_01 <- read_csv("~/Desktop/Postdoc/changepoint detection/IGRA/IGRA_1991_01.csv", col_names = FALSE)
dat <- IGRA_1991_01 %>% filter(X3 == as.Date('1991-01-10'), X1 == 7.3333, X2 == 134.4833, X4 == 11) %>% 
  arrange(X8) %>%
  dplyr::select(X6, X7, X8)
dat$X6 <- as.numeric(dat$X6)
dat$X7 <- as.numeric(dat$X7)
dat$X8 <- as.numeric(dat$X8)
colnames(dat) <- c('pressure', 'temperature','time')
xgrid = seq(max(log(dat$pressure)), min(log(dat$pressure)), length.out = length(log(dat$pressure)))
pred = approx(log(dat$pressure), dat$temperature, xout=xgrid)$y

# Create the plot
figure1a <- ggplot(dat, aes(x = log(pressure), y = temperature)) +
  # Line plot for raw data
  geom_line(aes(color = "Raw data"), size = 1.8) +
  # Points for raw data
  geom_point(aes(color = "Raw data"), size = 3) +
  # Line for predicted data (red, thinner)
  geom_line(data = data.frame(xgrid, pred), 
            aes(x = xgrid, y = pred, color = "Linear interpolation"), 
            size = 1.8) +  # Reduced line width
  # Vertical solid lines
  geom_vline(xintercept = log(200), linetype = "solid") +
  geom_vline(xintercept = log(1), linetype = "solid") +
  # Add text annotations
  annotate("text", x = 6.3, y = 315, label = "Troposphere", size = 5) +
  annotate("text", x = 3, y = 315, label = "Stratosphere", size = 5) +
  # Customize Legend
  scale_color_manual(
    name = "Legend",
    values = c("Raw data" = "black", "Linear interpolation" = "red")
  ) +
  # Axis and theme adjustments
  scale_x_reverse(limits = c(7, 0)) +  # Reverse x-axis to match xlim
  scale_y_continuous(limits = c(min(dat$temperature), 320)) +
  labs(x = "Log pressure (hPa)", y = "Temperature (K)") +  # Updated axis labels
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),  # Axis labels font size
    axis.text = element_text(size = 14, ),  # Tick mark font size
    axis.ticks = element_line(size = 1),  # Increase tick mark thickness
    axis.ticks.length = unit(0.25, "cm"),  # Increase tick length
    legend.position = c(0.6, 0.6),  # Legend inside plot
    legend.background = element_rect(fill = "transparent", color = NA),  # No border
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Legend text font size
    panel.grid = element_blank(),  # Remove all grid lines
    panel.border = element_blank(),  # Remove panel border
    axis.line = element_line(size = 1)  # Add x and y axis lines
  )

save_plot("Figure1a_new.png", figure1a, base_height = 3, base_width = 5,bg = "white")

# (b)

load("~/Desktop/Postdoc/changepoint detection/Main - Functional Changepoint Detection Code-selected/IGRA_global/IGRA_global_aggregate1/data.RData")
unique_latlon <- data %>% group_by(lat, long) %>% count()
dis_to_mt <- unique_latlon %>% 
  mutate(dis = sqrt((lat - 15.1429)^2 + (long - 120.3496)^2)) %>% 
  arrange(dis)
library(IndexNumR)
dat <- data %>% filter(lat == 7.3333, long == 134.4833)
dat$date <- as.Date(dat$date)
dat$pressure_level <- as.numeric(dat$pressure_level)
dat$temperature <- as.numeric(dat$temperature)
dat <- dat %>% dplyr::select(date, pressure_level, temperature, date_time)
dat <- dat[complete.cases(dat),]
dat_week <- dat %>% 
  mutate(week_index = weekIndex(as.Date(dat$date_time,format = "%Y-%m-%d"))) %>%
  group_by(week_index) %>% 
  mutate(temp = temperature, pressure = pressure_level)

unique(dat_week$week_index)
max(dat_week$pressure)
min(dat_week$pressure)
library(doParallel)
library(Hmisc)
dat_inte <- do.call(rbind, mclapply(1:53, function(xx){
  dat <- dat_week %>% filter(week_index == xx) %>% 
    group_by(pressure) %>% 
    summarise(week_index = xx, temp_unique = mean(temp, na.rm = TRUE)) %>% 
    arrange(pressure)
  smoothingSpline = smooth.spline(dat$pressure, dat$temp_unique, spar=0.35)
  temp_inte <- predict(smoothingSpline, c(seq(3.4,200,length.out = 5e4),
                                          seq(200, 1011.2, length.out = 1e4)))$y
  data.frame(week_index = xx, 
             pressure = c(seq(3.4,200,length.out = 5e4),
                          seq(200, 1011.2, length.out = 1e4)), 
             temp = temp_inte)
}, mc.cores = 7))

library(MASS)
library(scales)
rev_sqrt_trans <- function() {
  scales::trans_new(
    name = "rev_sqrt", 
    transform = function(x) -log10(x), 
    inverse = function(x) 10^x);
}
library(latex2exp)
figure1b <- ggplot(dat_inte, aes(x = week_index, y = log(pressure), fill = temp)) +
  # Heatmap tiles
  geom_tile(height = 0.1) + 
  
  # Color gradient for temperature
  scale_fill_gradientn(
    colors = c("blue", "lightblue", "yellow", "orange", "red"), 
    name = 'K'
    #name = TeX(r'(K($^{\circ}$F))')
  ) +
  
  # Reverse y-axis and define custom labels with no padding
  scale_y_continuous(
    trans = "reverse", 
    labels = c('6', '5', '4', '3', '2', '1.22'), 
    breaks = c(6, 5, 4, 3, 2, 1.22),
    expand = c(0, 0)  # Remove padding on the y-axis
  ) +
  
  # Custom x-axis labels with no padding
  scale_x_continuous(
    labels = c('0', '10', '20', '30', '40', '50'), 
    breaks = c(0, 10, 20, 30, 40, 50),
    expand = c(0, 0)  # Remove padding on the x-axis
  ) +
  
  # Vertical reference line
  geom_vline(xintercept = 24, linetype = "dashed", color = "black", size = 0.5) +
  
  # Axis labels
  labs(x = "Week", y = "Log pressure (hPa)") +
  
  # Theme customization
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),  # Axis label font size
    axis.text = element_text(size = 14),  # Tick mark font size
    axis.ticks = element_line(size = 1),  # Tick mark thickness
    axis.ticks.length = unit(0.25, "cm"),  # Tick length
    legend.text = element_text(size = 14),  # Legend text size
    legend.title = element_text(size = 14),  # Legend title size
    legend.key.size = unit(0.2, "lines"),  # Legend key size
    legend.key.height = unit(0.6, "cm"),  # Legend height
    legend.key.width = unit(0.2, "cm"),  # Legend width
    legend.background = element_blank(),  # Remove legend border
    legend.position = 'right',  # Adjust as needed
    panel.grid = element_blank(),  # Remove grid lines
    panel.border = element_blank(),  # Remove panel border
    axis.line = element_line(size = 1),  # Add x and y axis lines
    legend.spacing = unit(0, "cm"),  # Reduce the space between legend and plot
    legend.margin = margin(0, 0, 0, 0)  # Remove extra margin around legend
  )

save_plot("Figure1b_new.png", figure1b, base_height = 3, base_width = 5,bg = "white")

# (c)
max_lat = 90; max_long = 180
setwd("~/Desktop/Postdoc/changepoint detection/Main - Functional Changepoint Detection Code-selected")
load("grids_65_10_120.RData")
setwd("~/Desktop/Postdoc/changepoint detection/Main - Functional Changepoint Detection Code-selected/IGRA_global/IGRA_global_aggregate11_cstd")
load("select_grids.RData")
library(tidyverse)
library(ggplot2)
library(cowplot)
world <- map_data("world")
res = grids[select_grids, c('id','lat_lower','lat_upper','long_lower','long_upper')]
blocklist = data.frame(long_lower=c(-180,-180,-180,60,60), 
                       long_upper=c(-60,-60,-60,180,180), 
                       lat_lower=c(-5,-55,-65,-55,-65), 
                       lat_upper=c(5,-45,-55,-45,-55))
figure1c <- ggplot() +
  coord_fixed(ratio = 1.4) +
  
  # Draw rectangles for regions
  geom_rect(data = res, aes(xmin = long_lower, xmax = long_upper, ymin = lat_lower, ymax = lat_upper),
            colour = "black", fill = 'slategray1', size = 0.2) +
  
  # Draw rectangles for blocked regions
  geom_rect(data = blocklist, aes(xmin = long_lower, xmax = long_upper, ymin = lat_lower, ymax = lat_upper),
            colour = "black", fill = 'grey', size = 0) +
  
  # Draw world map
  geom_map(data = world, map = world, aes(long, lat, map_id = region),
           color = "black", fill = NA, size = 0.1) +
  
  # Set axis limits
  ylim(c(-max_lat, max_lat)) + xlim(c(-max_long, max_long)) +
  
  # Add points for locations
  geom_point(data = locations, aes(long, lat), size = 0.5, col = 'orange') +
  
  # Add points for specific locations
  geom_point(aes(y = 15.1383, x = 120.3500), colour = "red", shape = 17, size = 2) +
  geom_point(aes(y = 7.3333, x = 134.4833), colour = "blue", shape = 8, size = 2) +
  
  # Add horizontal and vertical lines
  geom_hline(yintercept = c(-55, -65), size = 0.2) +
  geom_vline(xintercept = c(-60, 60), size = 0.2) +
  
  # Remove axis labels
  xlab("") + ylab("") +
  
  # Add labels to regions
  geom_text(data = res, aes(long_lower + 10, (lat_lower + lat_upper) / 2, label = id, fontface = 2),
            hjust = 0.6, vjust = 0.6, size = 5, col = 1) + 
  
  # Set x-axis scale with custom labels
  scale_x_continuous(
    limits = c(-180, 180),
    breaks = c(-180, -60, 60, 180),
    labels = c("180\u00B0W", "60\u00B0W", "60\u00B0E", "180\u00B0E"),
    expand = c(0, 0),
    name = "",
    sec.axis = dup_axis()
  ) +
  
  # Set y-axis scale with custom labels
  scale_y_continuous(
    breaks = c(seq(-65, -5, 20), seq(5, 65, 20)),
    labels = c("65\u00B0S", "45\u00B0S", "25\u00B0S", "5\u00B0S",
               "5\u00B0N", "25\u00B0N", "45\u00B0N", "65\u00B0N"),
    expand = c(0, 0),
    name = "",
    sec.axis = dup_axis()
  ) +
  
  # Customize the theme
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.background = element_blank(),
    axis.text = element_text(size = 14),  # Same size as previous figure
    axis.title.x = element_blank(),  # Remove x-axis label
    axis.title.y = element_blank(),  # Remove y-axis label
    axis.ticks = element_line(size = 1),  # Same size as previous figure
    axis.ticks.length = unit(0.25, "cm"),  # Same length as previous figure
    axis.text.x = element_text(size = 14), # Bottom X-axis coordinates only
    axis.text.y = element_text(size = 14), # Right Y-axis coordinates only
    axis.text.x.top = element_blank(),    # Remove coordinates from the top
    axis.text.y.left = element_blank()    # Remove coordinates from the left
  )
setwd("~/Desktop/Postdoc/changepoint detection/modifacation")
save_plot("Figure1c_new.png", figure1c, base_height = 3, base_width = 5,bg = "white")


left_column <- plot_grid(
  figure1a, figure1b,    # Vertically stacked
  ncol = 1,              # 1 column in this grid
  rel_heights = c(1, 1.2),  # Adjust the height of each plot
  rel_widths = c(1, 1.5), # Adjust the relative widths
  labels = c("(a)", "(b)"),    # Optional labels for each plot
  label_x = c(0.07, 0.07),    # Position labels on the x-axis (first two labels closer to the left, last one on the right)
  label_y = c(1, 1.1)
)

# Combine the left_column grid with figure1c on the right
combined_plot <- plot_grid(
  left_column, figure1c, # Left column and the right plot
  ncol = 2,              # 2 columns in the main grid
  rel_widths = c(1, 1.1), # Adjust the relative widths
  rel_heights = c(1, 2),  # Adjust the height of each plot
  labels = c("", "(c)")    # Label for the right plot
)

figure1a_adjusted <- ggdraw(figure1a) + draw_label("(a)", x = 0.1, y = 1)
figure1b_adjusted <- ggdraw(figure1b) + draw_label("(b)", x = 0.1, y = 1.05)
figure1c_adjusted <- ggdraw(figure1c) + draw_label("(c)", x = 0, y = 0.75)

left_column <- plot_grid(
  ggdraw() + draw_plot(figure1a_adjusted, x = 0, y = 0, width = 0.89, height = 0.8),
  ggdraw() + draw_plot(figure1b_adjusted, x = 0, y = 0.1, width = 1, height = 0.9),
  nrow = 2,
  rel_heights = c(1, 1.2),
  align = "v",              # Align vertically
  axis = "l"                # Align left
)

# Combine the left_column grid with figure1c on the right
combined_plot <- plot_grid(
  left_column, 
  ggdraw() + draw_plot(figure1c_adjusted, x = 0, y = -0.25, width = 1, height = 1.5),
  ncol = 2,              # 2 columns in the main grid
  rel_widths = c(1, 1.1), # Adjust the relative widths
  rel_heights = c(1, 1),  # Adjust the height of each plot
  labels = c("", "")    # Label for the right plot
)

setwd("~/Desktop/Postdoc/changepoint detection/modifacation")
save_plot("Figure1_new.png", combined_plot, base_height = 6, base_width = 12,bg = "white")
##############

##############
## Figure 2
# (a)
load("~/Desktop/Postdoc/changepoint detection/piecewise/S_sim.RData")
load("~/Desktop/Postdoc/changepoint detection/piecewise/k_star.RData")
setwd("~/Desktop/Postdoc/changepoint detection/piecewise")
alpha.mean <- read.csv('alpha.mean_quadratic.csv')[,-1]
tau.mean <- read.csv('tau.mean_quadratic.csv')[,-1]
p.mean <- read.csv('p.mean_quadratic.csv')[,-1]
c.mean <- read.csv('c.mean_quadratic.csv')[,-1]
iters = 20000; thin_step = 10
r_thin = seq(iters*3/4+1, iters, thin_step) 
ns = 50
i = 30
N = nt = 50; x = (0:N)/N; nw = N+1
x <- seq(0,1,0.02)
nt = 50; nsim = 5; D = 21; x0 = (0: nt)/nt; D0 = (D-1)/2; Ds = rep(1:D0, each = 2)
sigma_delay = c(1, 1/Ds^3)
library(latex2exp)
mean_try = rep(alpha.mean, nw)*x*(1-x)+(ifelse(0<rep(p.mean[i],nw),1,0))*
  (ifelse(x<=rep(c.mean[i],nw),1,0))*N*x^2*rep(tau.mean[i],nw)*(1-rep(c.mean[i],nw))^2+
  (ifelse(0<rep(p.mean[i],nw),1,0))*(ifelse(x>rep(c.mean[i],nw),1,0))*N*
  (1-x)^2*rep(tau.mean[i],nw)*rep(c.mean[i],nw)^2
possible_change = seq(0, 1, 0.01)
possible_beta = seq(-1, -40, by = -0.01)
result_linear = matrix(0, ncol = 3, nrow = length(possible_beta)*length(possible_change))
colnames(result_linear) = c("squared error", "cp", "slope")
index = 1
for(ci in possible_change){
  for(betai in possible_beta){
    result_linear[index, 2:3] = c(ci, betai)
    mean = betai*((ci-1)*x0+(x0-ci)*(ifelse(x0>ci,1,0)))
    result_linear[index, 1] = sum((S_sim[, i] - mean)^2)
    index = index + 1
  }
}
cp_linear = result_linear[which.min(result_linear[, 1]), ][2]
beta_linear = result_linear[which.min(result_linear[, 1]), ][3]
x <- seq(0,1,0.02)
m_red = beta_linear*((cp_linear-1)*x+(x-cp_linear)*(ifelse(x>cp_linear,1,0)))

x <- seq(0, 1, 0.02)
time <- 0:N
i = 30

plot_data <- data.frame(
  Time = time,
  Y_TK = S_sim[, i],
  Mean_Try = mean_try,
  M_Red = m_red
)

plot_data <- data.frame(
  time = c(rep(time, 3),rep(k_star[i] * N,51)),
  dat = c(S_sim[, i], mean_try, m_red, seq(min(S_sim[, i]), max(S_sim[, i]), length.out = 51)),
  type = rep(c('Y', 'quadratic', 'linear', 'true'), each = 51)
)

figure2a <- ggplot(plot_data, aes(x=time, y=dat)) +
  geom_line(aes(color=type), size = 1)+
  geom_vline(xintercept = c.mean[i] * N, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = cp_linear * N, color = "green", linetype = "dashed", size = 1) +
  # Custom Colors and Legend
  scale_color_manual(
    values = c('green','red','black', 'blue'),
    name = NULL,
    labels = c("Linear Fit", 
               "Quadratic Fit", 
               "True Change Point", 
               TeX(r'($Y_{T,k}$)'))
  ) +
  scale_linetype_identity(guide = "none") +
  # Axis Labels
  labs(
    x = "Time",
    y = TeX(r'($Y_{T,k}$ process)')
  ) +
  # Theme for Aesthetic Adjustments
  theme_classic() +  # Removes gridlines, keeps axes visible
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.ticks = element_line(size = 1),  # Tick mark thickness
    axis.ticks.length = unit(0.25, "cm"),  # Tick length
    legend.text = element_text(size = 14),
    legend.position = c(0.75, 0.85),  # Position legend inside the plot
    legend.background = element_blank(),  # Remove legend background
    legend.key = element_blank(),         # Remove legend key background
    legend.margin = margin(t = 4, r = 4, b = 4, l = 4)
  )

setwd("~/Desktop/Postdoc/changepoint detection/modifacation")
save_plot("Figure2a_new.png", figure2a, base_height = 3, base_width = 5,bg = "white")

# (b)
tau_index <- data.frame(tau.mean = tau.mean, 
                        index = c(1:50))
tau_index <- tau_index %>% mutate(null_alt = ifelse(index %in% c(2,3,7,8,38), 'Null','Alternative'))
figure2b <- ggplot(tau_index, aes(x = tau.mean, color = null_alt, fill = null_alt)) +
  geom_histogram(position = "identity", bins = 30, color = "black") +
  scale_fill_manual(values = c("Null" = "darkorange", "Alternative" = "deepskyblue")) +
  scale_color_manual(values = c("Null" = "darkorange", "Alternative" = "deepskyblue")) +
  labs(
    x = expression(tau(bold(s))),
    y = "Frequency",
    fill = NULL,  # Removes legend title
    color = NULL  # Removes legend title
  ) +
  theme_classic() +  # Removes grid background
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.ticks = element_line(size = 1),  # Tick mark thickness
    axis.ticks.length = unit(0.25, "cm"),  # Tick length
    legend.text = element_text(size = 12),
    legend.position = c(0.2, 0.8),  # Adjust legend position (x, y)
    legend.background = element_blank(),  # Removes legend background
    legend.key = element_rect(fill = "transparent", color = "transparent")
  )

setwd("~/Desktop/Postdoc/changepoint detection/modifacation")
save_plot("Figure2b_new.png", figure2b, base_height = 3, base_width = 5,bg = "white")

figure2a_adjusted <- ggdraw(figure2a) + draw_label("(a)", x = 0.05, y = 0.95)
figure2b_adjusted <- ggdraw(figure2b) + draw_label("(b)", x = 0.02, y = 0.982)

combined_plot <- plot_grid(
  ggdraw() + draw_plot(figure2a_adjusted, x = 0, y = 0),
  ggdraw() + draw_plot(figure2b_adjusted, x = 0, y = -0.032),
  ncol = 2,              # 2 columns in the main grid
  labels = c("", "")     # Label for the right plot
)
setwd("~/Desktop/Postdoc/changepoint detection/modifacation")
save_plot("Figure2_new.png", combined_plot, base_height = 3, base_width = 10,bg = "white")
##############

##############
## Figure 3

library(truncnorm)
set.seed(1234)
original_dist <- rtruncnorm(200, a=0, b=Inf, mean = 0, sd = 1)
df <- data.frame(
  Data=factor(rep(c("g=1", "g=10"), each=200)),
  weight=c(original_dist, 
           rtruncnorm(200, a=0, b=Inf, mean = 0, sd = 10))
)

d1dens <- with(df, density(weight[Data == "g=1"], 
                           from = min(weight), 
                           to = max(weight)))
d2dens <- with(df, density(weight[Data == "g=10"], 
                           from = min(weight),
                           to = max(weight)))
joint <- pmin(d1dens$y, d2dens$y)
df2 <- data.frame(x = rep(d1dens$x, 3), 
                  y = c(d1dens$y, d2dens$y, joint),
                  Distribution = rep(c("original distribution", "g=10", "Overlap"), each = length(d1dens$x)))
library(ggplot2)
g10 = ggplot(df2, aes(x, y, fill = Distribution)) + 
  scale_fill_discrete(labels=c("Slab","Spike", "Overlap"))+
  scale_y_sqrt() + 
  geom_area(position = position_identity(), color = "black", stat = 'identity')+
  theme_bw() +
  annotate(geom="text", x=15, y=0.75, label= 'g=10', size = 5)+
  theme_classic() +  # Removes grid background
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.ticks = element_line(size = 1),  # Tick mark thickness
    axis.ticks.length = unit(0.25, "cm"),  # Tick length
    legend.text = element_text(size = 14),
    legend.position = c(0.8, 0.8),  # Adjust legend position (x, y)
    legend.background = element_rect(fill = "white"),  # Removes legend background
    legend.key = element_rect(fill = "transparent", color = "transparent")
  )
#annotate(geom="text", x=15, y=0.75, label= paste('overlap = ',round(sum(joint) / sum(d1dens$y, d2dens$y), 5)))

df <- data.frame(
  Data=factor(rep(c("g=1", "g=100"), each=200)),
  weight=c(original_dist, 
           rtruncnorm(200, a=0, b=Inf, mean = 0, sd = 100))
)

d1dens <- with(df, density(weight[Data == "g=1"], 
                           from = min(weight), 
                           to = max(weight)))
d2dens <- with(df, density(weight[Data == "g=100"], 
                           from = min(weight),
                           to = max(weight)))
joint <- pmin(d1dens$y, d2dens$y)
df2 <- data.frame(x = rep(d1dens$x, 3), 
                  y = c(d1dens$y, d2dens$y, joint),
                  Distribution = rep(c("original distribution", "g=100", "Overlap"), each = length(d1dens$x)))
library(ggplot2)
g100 <- ggplot(df2, aes(x, y, fill = Distribution)) + 
  scale_fill_discrete(labels=c("Slab", "Spike", "Overlap")) +
  scale_y_sqrt()+ 
  geom_area(position = position_identity(), color = "black") +
  theme_bw() + 
  xlim(c(0,50)) +
  annotate(geom="text", x=25, y=0.75, label= 'g=100', size = 5)+
  theme_classic() +  # Removes grid background
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.ticks = element_line(size = 1),  # Tick mark thickness
    axis.ticks.length = unit(0.25, "cm"),  # Tick length
    legend.text = element_text(size = 14),
    legend.position = c(0.8, 0.8),  # Adjust legend position (x, y)
    legend.background = element_rect(fill = "white"),  # Removes legend background
    legend.key = element_rect(fill = "transparent", color = "transparent")
  )
# annotate(geom="text", x=25, y=0.75, label= paste('overlap = ',round(sum(joint) / sum(d1dens$y, d2dens$y), 5)))

df <- data.frame(
  Data=factor(rep(c("g=1", "g=1000"), each=200)),
  weight=c(original_dist, 
           rtruncnorm(200, a=0, b=Inf, mean = 0, sd = 1000))
)

d1dens <- with(df, density(weight[Data == "g=1"], 
                           from = min(weight), 
                           to = max(weight)))
d2dens <- with(df, density(weight[Data == "g=1000"], 
                           from = min(weight),
                           to = max(weight)))
joint <- pmin(d1dens$y, d2dens$y)
df2 <- data.frame(x = rep(d2dens$x, 3), 
                  y = c(d1dens$y, d2dens$y, joint),
                  Distribution = rep(c("original distribution", "g=1000", "Overlap"), each = length(d1dens$x)))
library(ggplot2)
g1000 = ggplot(df2, aes(x, y, fill = Distribution)) + 
  scale_fill_discrete(labels=c("Slab", "Spike", "Overlap")) +
  scale_y_sqrt() + 
  geom_area(position = position_identity(), color = "black") +
  theme_bw() +
  xlim(c(0,100)) + 
  annotate(geom="text", x=50, y=1.3, label= 'g=1000', size = 5)+
  theme_classic() +  # Removes grid background
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.ticks = element_line(size = 1),  # Tick mark thickness
    axis.ticks.length = unit(0.25, "cm"),  # Tick length
    legend.text = element_text(size = 14),
    legend.position = c(0.8, 0.8),  # Adjust legend position (x, y)
    legend.background = element_rect(fill = "white"),  # Removes legend background
    legend.key = element_rect(fill = "transparent", color = "transparent")
  )
# annotate(geom="text", x=50, y=1.3, label= paste('overlap = ',round(sum(joint) / sum(d1dens$y, d2dens$y), 5)))

library(ggpubr)
figure3 <- ggarrange(g10, g100, g1000, 
          labels = c("(a)", "(b)", "(c)"),
          hjust=0.05,
          ncol = 3, nrow = 1, 
          common.legend = TRUE, legend="bottom")
setwd("~/Desktop/Postdoc/changepoint detection/modifacation")
save_plot("Figure3_new.png", figure3, base_height = 4, base_width = 10,bg = "white")
##############

##############
## Figure 6
library(scales)
setwd("~/Desktop/Postdoc/changepoint detection/geosphere")
Nt = 43
prec_temp = read.csv("qnormc.csv", header = TRUE)[,-1]
r_thin = seq(15001, 20000, 10)
c.mean = apply(pnorm(as.matrix(prec_temp[r_thin,])), 2, mean)
xgrid = seq(300, 20, -10)
cp_BH = round(as.vector(c.mean)*43)
load("~/Desktop/Postdoc/changepoint detection/Main - Functional Changepoint Detection Code-selected/IGRA_global/IGRA_global_aggregate11_cstd/select_grids.RData")
load("~/Desktop/Postdoc/changepoint detection/Main - Functional Changepoint Detection Code-selected/IGRA_global/IGRA_global_aggregate11_cstd/cstd_data_list.RData")

i = which(select_grids == 8)
# Sample Data (Replace with your actual data)
obs <- cstd_data_list[[i]]
cp_BH_i <- cp_BH[i]
select_grid <- select_grids[i]

# Convert to long format for ggplot
library(tidyr)
obs_long <- data.frame(xgrid = rep(xgrid, Nt), 
                       temperature = as.vector(obs), 
                       time = rep(1:Nt, each = nrow(obs)))

obs_long$color_group <- ifelse(obs_long$time <= cp_BH_i, "Before CP", "After CP")

# Calculate mean before and after change point

mean_data_before <- data.frame(
  xgrid = xgrid,
  temperature = rowMeans(obs[, 1:cp_BH[i]]),  # Assuming cp_BH[i] marks the breakpoint
  time = Nt+1
)
mean_data_after <- data.frame(
  xgrid = xgrid,
  temperature = rowMeans(obs[, (cp_BH_i+1):Nt]),  # Assuming cp_BH[i] marks the breakpoint
  time = Nt+2
)

# Plot
figure6a1 <- ggplot(obs_long, aes(x = xgrid, y = temperature, group = time)) +
  # Plot individual lines
  geom_line(aes(color = color_group)) +
  
  # Add mean line before change point
  geom_line(data = mean_data_before, aes(x = xgrid, y = temperature), 
            color = 'blue', size = 1) +
  
  # Add mean line after change point
  geom_line(data = mean_data_after, aes(x = xgrid, y = temperature), 
            color = 'red', size = 1) +
  
  # Manual colors for groups
  scale_color_manual(values = c("Before CP" = alpha("blue", 0.35), 
                                "After CP" = alpha("red", 0.35))) +
  
  # Reverse x-axis
  scale_x_reverse() +
  
  # Labels
  labs(
    x = "Pressure level (hPa)", 
    y = "Temperature Anomaly"
  ) +
  
  # Clean theme
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.ticks = element_line(size = 1),  # Tick mark thickness
    axis.ticks.length = unit(0.25, "cm"),  # Tick length
    
    # Remove legend
    legend.position = "none"
  ) +
  
  # Annotate grid label
  annotate("text", x = max(xgrid) - 25, y = max(obs) + 1, 
           label = paste('Grid ', select_grid,' (10° N, 120° W)', sep = ''), 
           hjust = 0, size = 5)

i = which(select_grids == 20)
# Sample Data (Replace with your actual data)
obs <- cstd_data_list[[i]]
cp_BH_i <- cp_BH[i]
select_grid <- select_grids[i]

# Convert to long format for ggplot
library(tidyr)
obs_long <- data.frame(xgrid = rep(xgrid, Nt), 
                       temperature = as.vector(obs), 
                       time = rep(1:Nt, each = nrow(obs)))

obs_long$color_group <- ifelse(obs_long$time <= cp_BH_i, "Before CP", "After CP")

# Calculate mean before and after change point

mean_data_before <- data.frame(
  xgrid = xgrid,
  temperature = rowMeans(obs[, 1:cp_BH[i]]),  # Assuming cp_BH[i] marks the breakpoint
  time = Nt+1
)
mean_data_after <- data.frame(
  xgrid = xgrid,
  temperature = rowMeans(obs[, (cp_BH_i+1):Nt]),  # Assuming cp_BH[i] marks the breakpoint
  time = Nt+2
)

# Plot
figure6a2 <- ggplot(obs_long, aes(x = xgrid, y = temperature, group = time)) +
  # Plot individual lines
  geom_line(aes(color = color_group)) +
  
  # Add mean line before change point
  geom_line(data = mean_data_before, aes(x = xgrid, y = temperature), 
            color = 'blue', size = 1) +
  
  # Add mean line after change point
  geom_line(data = mean_data_after, aes(x = xgrid, y = temperature), 
            color = 'red', size = 1) +
  
  # Manual colors for groups
  scale_color_manual(values = c("Before CP" = alpha("blue", 0.35), 
                                "After CP" = alpha("red", 0.35))) +
  
  # Reverse x-axis
  scale_x_reverse() +
  
  # Labels
  labs(
    x = "Pressure level (hPa)", 
    y = "Temperature Anomaly"
  ) +
  
  # Clean theme
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.ticks = element_line(size = 1),  # Tick mark thickness
    axis.ticks.length = unit(0.25, "cm"),  # Tick length
    
    # Remove legend
    legend.position = "none"
  ) +
  
  # Annotate grid label
  annotate("text", x = max(xgrid) - 25, y = max(obs) + 1, 
           label = paste('Grid ',select_grid,' (0°, 0° )', sep = ''), 
           hjust = 0, size = 5)

i = which(select_grids == 32)
# Sample Data (Replace with your actual data)
obs <- cstd_data_list[[i]]
cp_BH_i <- cp_BH[i]
select_grid <- select_grids[i]

# Convert to long format for ggplot
library(tidyr)
obs_long <- data.frame(xgrid = rep(xgrid, Nt), 
                       temperature = as.vector(obs), 
                       time = rep(1:Nt, each = nrow(obs)))

obs_long$color_group <- ifelse(obs_long$time <= cp_BH_i, "Before CP", "After CP")

# Calculate mean before and after change point

mean_data_before <- data.frame(
  xgrid = xgrid,
  temperature = rowMeans(obs[, 1:cp_BH[i]]),  # Assuming cp_BH[i] marks the breakpoint
  time = Nt+1
)
mean_data_after <- data.frame(
  xgrid = xgrid,
  temperature = rowMeans(obs[, (cp_BH_i+1):Nt]),  # Assuming cp_BH[i] marks the breakpoint
  time = Nt+2
)

# Plot
figure6a3 <- ggplot(obs_long, aes(x = xgrid, y = temperature, group = time)) +
  # Plot individual lines
  geom_line(aes(color = color_group)) +
  
  # Add mean line before change point
  geom_line(data = mean_data_before, aes(x = xgrid, y = temperature), 
            color = 'blue', size = 1) +
  
  # Add mean line after change point
  geom_line(data = mean_data_after, aes(x = xgrid, y = temperature), 
            color = 'red', size = 1) +
  
  # Manual colors for groups
  scale_color_manual(values = c("Before CP" = alpha("blue", 0.35), 
                                "After CP" = alpha("red", 0.35))) +
  
  # Reverse x-axis
  scale_x_reverse() +
  
  # Labels
  labs(
    x = "Pressure level (hPa)", 
    y = "Temperature Anomaly"
  ) +
  
  # Clean theme
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.ticks = element_line(size = 1),  # Tick mark thickness
    axis.ticks.length = unit(0.25, "cm"),  # Tick length
    
    # Remove legend
    legend.position = "none"
  ) +
  
  # Annotate grid label
  annotate("text", x = max(xgrid) - 25, y = max(obs) + 1, 
           label = paste('Grid ',select_grid,' (10° S, 120° E)', sep = ''), 
           hjust = 0, size = 5)

figure6a <- plot_grid(figure6a1, figure6a2, figure6a3, ncol = 3)
setwd("~/Desktop/Postdoc/changepoint detection/modifacation")
save_plot("Figure6a_new.png", figure6a, base_height = 4, base_width = 10,bg = "white")

i = which(select_grids == 10)
# Sample Data (Replace with your actual data)
obs <- cstd_data_list[[i]]
cp_BH_i <- cp_BH[i]
select_grid <- select_grids[i]

# Convert to long format for ggplot
library(tidyr)
obs_long <- data.frame(xgrid = rep(xgrid, Nt), 
                       temperature = as.vector(obs), 
                       time = rep(1:Nt, each = nrow(obs)))

obs_long$color_group <- ifelse(obs_long$time <= cp_BH_i, "Before CP", "After CP")

# Calculate mean before and after change point

mean_data_before <- data.frame(
  xgrid = xgrid,
  temperature = rowMeans(obs[, 1:cp_BH[i]]),  # Assuming cp_BH[i] marks the breakpoint
  time = Nt+1
)
mean_data_after <- data.frame(
  xgrid = xgrid,
  temperature = rowMeans(obs[, (cp_BH_i+1):Nt]),  # Assuming cp_BH[i] marks the breakpoint
  time = Nt+2
)

# Plot
figure6b1 <- ggplot(obs_long, aes(x = xgrid, y = temperature, group = time)) +
  # Plot individual lines
  geom_line(aes(color = color_group)) +
  
  # Add mean line before change point
  geom_line(data = mean_data_before, aes(x = xgrid, y = temperature), 
            color = 'blue', size = 1) +
  
  # Add mean line after change point
  geom_line(data = mean_data_after, aes(x = xgrid, y = temperature), 
            color = 'red', size = 1) +
  
  # Manual colors for groups
  scale_color_manual(values = c("Before CP" = alpha("blue", 0.35), 
                                "After CP" = alpha("red", 0.35))) +
  
  # Reverse x-axis
  scale_x_reverse() +
  
  # Labels
  labs(
    x = "Pressure level (hPa)", 
    y = "Temperature Anomaly"
  ) +
  
  # Clean theme
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.ticks = element_line(size = 1),  # Tick mark thickness
    axis.ticks.length = unit(0.25, "cm"),  # Tick length
    
    # Remove legend
    legend.position = "none"
  ) +
  
  # Annotate grid label
  annotate("text", x = max(xgrid) - 25, y = max(obs) + 1, 
           label = paste('Grid ',select_grid,' (30° N, 120° W)', sep = ''), 
           hjust = 0, size = 5)

i = which(select_grids == 23)
# Sample Data (Replace with your actual data)
obs <- cstd_data_list[[i]]
cp_BH_i <- cp_BH[i]
select_grid <- select_grids[i]

# Convert to long format for ggplot
library(tidyr)
obs_long <- data.frame(xgrid = rep(xgrid, Nt), 
                       temperature = as.vector(obs), 
                       time = rep(1:Nt, each = nrow(obs)))

obs_long$color_group <- ifelse(obs_long$time <= cp_BH_i, "Before CP", "After CP")

# Calculate mean before and after change point

mean_data_before <- data.frame(
  xgrid = xgrid,
  temperature = rowMeans(obs[, 1:cp_BH[i]]),  # Assuming cp_BH[i] marks the breakpoint
  time = Nt+1
)
mean_data_after <- data.frame(
  xgrid = xgrid,
  temperature = rowMeans(obs[, (cp_BH_i+1):Nt]),  # Assuming cp_BH[i] marks the breakpoint
  time = Nt+2
)

# Plot
figure6b2 <- ggplot(obs_long, aes(x = xgrid, y = temperature, group = time)) +
  # Plot individual lines
  geom_line(aes(color = color_group)) +
  
  # Add mean line before change point
  geom_line(data = mean_data_before, aes(x = xgrid, y = temperature), 
            color = 'blue', size = 1) +
  
  # Add mean line after change point
  geom_line(data = mean_data_after, aes(x = xgrid, y = temperature), 
            color = 'red', size = 1) +
  
  # Manual colors for groups
  scale_color_manual(values = c("Before CP" = alpha("blue", 0.35), 
                                "After CP" = alpha("red", 0.35))) +
  
  # Reverse x-axis
  scale_x_reverse() +
  
  # Labels
  labs(
    x = "Pressure level (hPa)", 
    y = "Temperature Anomaly"
  ) +
  
  # Clean theme
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.ticks = element_line(size = 1),  # Tick mark thickness
    axis.ticks.length = unit(0.25, "cm"),  # Tick length
    
    # Remove legend
    legend.position = "none"
  ) +
  
  # Annotate grid label
  annotate("text", x = max(xgrid) - 25, y = max(obs) + 1, 
           label = paste('Grid ',select_grid,' (30° N, 0° )', sep = ''), 
           hjust = 0, size = 5)

i = which(select_grids == 36)
# Sample Data (Replace with your actual data)
obs <- cstd_data_list[[i]]
cp_BH_i <- cp_BH[i]
select_grid <- select_grids[i]

# Convert to long format for ggplot
library(tidyr)
obs_long <- data.frame(xgrid = rep(xgrid, Nt), 
                       temperature = as.vector(obs), 
                       time = rep(1:Nt, each = nrow(obs)))

obs_long$color_group <- ifelse(obs_long$time <= cp_BH_i, "Before CP", "After CP")

# Calculate mean before and after change point

mean_data_before <- data.frame(
  xgrid = xgrid,
  temperature = rowMeans(obs[, 1:cp_BH[i]]),  # Assuming cp_BH[i] marks the breakpoint
  time = Nt+1
)
mean_data_after <- data.frame(
  xgrid = xgrid,
  temperature = rowMeans(obs[, (cp_BH_i+1):Nt]),  # Assuming cp_BH[i] marks the breakpoint
  time = Nt+2
)

# Plot
figure6b3 <- ggplot(obs_long, aes(x = xgrid, y = temperature, group = time)) +
  # Plot individual lines
  geom_line(aes(color = color_group)) +
  
  # Add mean line before change point
  geom_line(data = mean_data_before, aes(x = xgrid, y = temperature), 
            color = 'blue', size = 1) +
  
  # Add mean line after change point
  geom_line(data = mean_data_after, aes(x = xgrid, y = temperature), 
            color = 'red', size = 1) +
  
  # Manual colors for groups
  scale_color_manual(values = c("Before CP" = alpha("blue", 0.35), 
                                "After CP" = alpha("red", 0.35))) +
  
  # Reverse x-axis
  scale_x_reverse() +
  
  # Labels
  labs(
    x = "Pressure level (hPa)", 
    y = "Temperature Anomaly"
  ) +
  
  # Clean theme
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.ticks = element_line(size = 1),  # Tick mark thickness
    axis.ticks.length = unit(0.25, "cm"),  # Tick length
    
    # Remove legend
    legend.position = "none"
  ) +
  
  # Annotate grid label
  annotate("text", x = max(xgrid) - 25, y = max(obs) + 1, 
           label = paste('Grid ',select_grid,' (30° N, 120° E)', sep = ''), 
           hjust = 0, size = 5)

figure6b <- plot_grid(figure6b1, figure6b2, figure6b3, ncol = 3)
setwd("~/Desktop/Postdoc/changepoint detection/modifacation")
save_plot("Figure6b_new.png", figure6b, base_height = 4, base_width = 10,bg = "white")

i = which(select_grids == 5)
# Sample Data (Replace with your actual data)
obs <- cstd_data_list[[i]]
cp_BH_i <- cp_BH[i]
select_grid <- select_grids[i]

# Convert to long format for ggplot
library(tidyr)
obs_long <- data.frame(xgrid = rep(xgrid, Nt), 
                       temperature = as.vector(obs), 
                       time = rep(1:Nt, each = nrow(obs)))

obs_long$color_group <- ifelse(obs_long$time <= cp_BH_i, "Before CP", "After CP")

# Calculate mean before and after change point

mean_data_before <- data.frame(
  xgrid = xgrid,
  temperature = rowMeans(obs[, 1:cp_BH[i]]),  # Assuming cp_BH[i] marks the breakpoint
  time = Nt+1
)
mean_data_after <- data.frame(
  xgrid = xgrid,
  temperature = rowMeans(obs[, (cp_BH_i+1):Nt]),  # Assuming cp_BH[i] marks the breakpoint
  time = Nt+2
)

# Plot
figure6c1 <- ggplot(obs_long, aes(x = xgrid, y = temperature, group = time)) +
  # Plot individual lines
  geom_line(aes(color = color_group)) +
  
  # Add mean line before change point
  geom_line(data = mean_data_before, aes(x = xgrid, y = temperature), 
            color = 'blue', size = 1) +
  
  # Add mean line after change point
  geom_line(data = mean_data_after, aes(x = xgrid, y = temperature), 
            color = 'red', size = 1) +
  
  # Manual colors for groups
  scale_color_manual(values = c("Before CP" = alpha("blue", 0.35), 
                                "After CP" = alpha("red", 0.35))) +
  
  # Reverse x-axis
  scale_x_reverse() +
  
  # Labels
  labs(
    x = "Pressure level (hPa)", 
    y = "Temperature Anomaly"
  ) +
  
  # Clean theme
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.ticks = element_line(size = 1),  # Tick mark thickness
    axis.ticks.length = unit(0.25, "cm"),  # Tick length
    
    # Remove legend
    legend.position = "none"
  ) +
  
  # Annotate grid label
  annotate("text", x = max(xgrid) - 25, y = max(obs) + 1, 
           label = paste('Grid ',select_grid,' (20° S, 120° W)', sep = ''), 
           hjust = 0, size = 5)

i = which(select_grids == 18)
# Sample Data (Replace with your actual data)
obs <- cstd_data_list[[i]]
cp_BH_i <- cp_BH[i]
select_grid <- select_grids[i]

# Convert to long format for ggplot
library(tidyr)
obs_long <- data.frame(xgrid = rep(xgrid, Nt), 
                       temperature = as.vector(obs), 
                       time = rep(1:Nt, each = nrow(obs)))

obs_long$color_group <- ifelse(obs_long$time <= cp_BH_i, "Before CP", "After CP")

# Calculate mean before and after change point

mean_data_before <- data.frame(
  xgrid = xgrid,
  temperature = rowMeans(obs[, 1:cp_BH[i]]),  # Assuming cp_BH[i] marks the breakpoint
  time = Nt+1
)
mean_data_after <- data.frame(
  xgrid = xgrid,
  temperature = rowMeans(obs[, (cp_BH_i+1):Nt]),  # Assuming cp_BH[i] marks the breakpoint
  time = Nt+2
)

# Plot
figure6c2 <- ggplot(obs_long, aes(x = xgrid, y = temperature, group = time)) +
  # Plot individual lines
  geom_line(aes(color = color_group)) +
  
  # Add mean line before change point
  geom_line(data = mean_data_before, aes(x = xgrid, y = temperature), 
            color = 'blue', size = 1) +
  
  # Add mean line after change point
  geom_line(data = mean_data_after, aes(x = xgrid, y = temperature), 
            color = 'red', size = 1) +
  
  # Manual colors for groups
  scale_color_manual(values = c("Before CP" = alpha("blue", 0.35), 
                                "After CP" = alpha("red", 0.35))) +
  
  # Reverse x-axis
  scale_x_reverse() +
  
  # Labels
  labs(
    x = "Pressure level (hPa)", 
    y = "Temperature Anomaly"
  ) +
  
  # Clean theme
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.ticks = element_line(size = 1),  # Tick mark thickness
    axis.ticks.length = unit(0.25, "cm"),  # Tick length
    
    # Remove legend
    legend.position = "none"
  ) +
  
  # Annotate grid label
  annotate("text", x = max(xgrid) - 25, y = max(obs) + 1, 
           label = paste('Grid ',select_grid,' (20° S, 0° )', sep = ''), 
           hjust = 0, size = 5)

i = which(select_grids == 31)
# Sample Data (Replace with your actual data)
obs <- cstd_data_list[[i]]
cp_BH_i <- cp_BH[i]
select_grid <- select_grids[i]

# Convert to long format for ggplot
library(tidyr)
obs_long <- data.frame(xgrid = rep(xgrid, Nt), 
                       temperature = as.vector(obs), 
                       time = rep(1:Nt, each = nrow(obs)))

obs_long$color_group <- ifelse(obs_long$time <= cp_BH_i, "Before CP", "After CP")

# Calculate mean before and after change point

mean_data_before <- data.frame(
  xgrid = xgrid,
  temperature = rowMeans(obs[, 1:cp_BH[i]]),  # Assuming cp_BH[i] marks the breakpoint
  time = Nt+1
)
mean_data_after <- data.frame(
  xgrid = xgrid,
  temperature = rowMeans(obs[, (cp_BH_i+1):Nt]),  # Assuming cp_BH[i] marks the breakpoint
  time = Nt+2
)

# Plot
figure6c3 <- ggplot(obs_long, aes(x = xgrid, y = temperature, group = time)) +
  # Plot individual lines
  geom_line(aes(color = color_group)) +
  
  # Add mean line before change point
  geom_line(data = mean_data_before, aes(x = xgrid, y = temperature), 
            color = 'blue', size = 1) +
  
  # Add mean line after change point
  geom_line(data = mean_data_after, aes(x = xgrid, y = temperature), 
            color = 'red', size = 1) +
  
  # Manual colors for groups
  scale_color_manual(values = c("Before CP" = alpha("blue", 0.35), 
                                "After CP" = alpha("red", 0.35))) +
  
  # Reverse x-axis
  scale_x_reverse() +
  
  # Labels
  labs(
    x = "Pressure level (hPa)", 
    y = "Temperature Anomaly"
  ) +
  
  # Clean theme
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.ticks = element_line(size = 1),  # Tick mark thickness
    axis.ticks.length = unit(0.25, "cm"),  # Tick length
    
    # Remove legend
    legend.position = "none"
  ) +
  
  # Annotate grid label
  annotate("text", x = max(xgrid) - 25, y = max(obs) + 1, 
           label = paste('Grid ',select_grid,' (20° S, 120° E)', sep = ''), 
           hjust = 0, size = 5)

figure6c <- plot_grid(figure6c1, figure6c2, figure6c3, ncol = 3)
setwd("~/Desktop/Postdoc/changepoint detection/modifacation")
save_plot("Figure6c_new.png", figure6c, base_height = 3, base_width = 9,bg = "white")

setwd("~/Desktop/Postdoc/changepoint detection/modifacation")
figure6 <- plot_grid(figure6a, figure6b, figure6c, nrow = 3, 
                     labels = c("(a)","(b)", "(c)"), label_size = 16)
save_plot("Figure6_new.png", figure6, base_height = 11, base_width = 13,bg = "white")

##############

##############
## Figure 5
#(a)
setwd("~/Desktop/Postdoc/changepoint detection/geosphere")
p_temp = read.csv("p.csv", header = TRUE)[,-1]
prec_temp = read.csv("qnormc.csv", header = TRUE)[,-1]
r_thin = seq(15001, 20000, 10)
c.mean = apply(pnorm(as.matrix(prec_temp[r_thin,])), 2, mean)
p.mean = apply(as.matrix(p_temp[r_thin,]), 2, mean)
max_lat = 65; max_long = 180
setwd("~/Desktop/Postdoc/changepoint detection/Main - Functional Changepoint Detection Code-selected")
load("grids_65_10_120.RData")
res = grids[select_grids, c('id','lat_lower','lat_upper','long_lower','long_upper')]
res$Changepoint = round(as.vector(c.mean)*43)
range(res$Changepoint)
br = seq(7,41,2)
lb = c("06/09/91", "06/23/91","07/07/91","07/21/91","08/04/91","08/18/91",
       "09/01/91", "09/15/91","09/29/91","10/13/91","10/27/91","11/10/91", 
       "11/24/91", "12/08/91","12/22/91","01/05/92","01/19/91","02/02/92")
res$Type = ifelse(as.vector(p.mean) > 0.4, "alt", "null")
res$Type = factor(res$Type, levels = c("alt", "null"))
res$altChangepoint = res$Changepoint
#res$altChangepoint[res$Type=="null"] = NA
range(res$altChangepoint, na.rm=T)
library(viridis)
hatch_data <- c()
for (i in 1:nrow(res)) {
  if (res$Type[i] == 'null') {
    hatch_data <- rbind(hatch_data, 
                        data.frame(x = rep(seq(res$long_lower[i], res$long_upper[i], by = 4), each = 1),
                                   xend = rep(seq(res$long_lower[i], res$long_upper[i], by = 4), each = 1),
                                   y = rep(res$lat_lower[i], length(seq(res$long_lower[i], res$long_upper[i], by = 4)) * 1),
                                   yend = rep(res$lat_upper[i], length(seq(res$long_lower[i], res$long_upper[i], by = 4)) * 1)))
  }
}
library(ggpattern)
figure5a = ggplot() +
  coord_fixed(ratio = 1.4) +
  geom_rect(data = res, aes(xmin = long_lower, xmax = long_upper, ymin = lat_lower, ymax = lat_upper, fill = altChangepoint)) +
  geom_map(data = world, map = world, aes(long, lat, map_id = region),
           color = "black", fill = NA, size = 0.1) +
  scale_fill_gradient2(breaks = br, limits = c(7, 41), labels = lb, high = "dodgerblue1", mid = "white", midpoint = 7,
                       na.value = "white", name = "Changepoint") + 
  geom_segment(data = hatch_data, aes(x = x, y = y, xend = xend, yend = yend),
               color = "grey15", size = 0.2, linetype = "solid") + 
  ylim(c(-max_lat, max_lat)) + xlim(c(-max_long, max_long)) +
  geom_point(aes(y = 15.1383, x = 120.3500), colour = "red", shape = 17, size = 2) +
  xlab("Longitude") + ylab("Latitude") +
  geom_text(data = grid_center[select_grids,], aes(long_lower + 10, (lat_lower + lat_upper) / 2, label = id, fontface = 2), 
            hjust = 0.6, vjust = 0.6, size = 2.5, col = 1) +
  theme(
    legend.text = element_text(size = 6), 
    legend.key.height = unit(1.1, 'cm'),
    legend.title = element_text(size = 8),
    axis.title = element_text(size = 8),
    axis.text.x.bottom = element_blank(),  # Remove bottom longitude labels
    axis.ticks.x.bottom = element_blank()  # Remove bottom longitude ticks
  ) + 
  scale_x_continuous(limits = c(-180, 180), breaks = c(-180, -60, 60, 180), 
                     labels = c("180\u00B0W", "60\u00B0W", "60\u00B0E", "180\u00B0E"), 
                     expand = c(0, 0), 
                     name = "",
                     sec.axis = dup_axis()) +
  scale_y_continuous(breaks = c(seq(-65, -5, 20), seq(5, 65, 20)), 
                     labels = c("65\u00B0S", "45\u00B0S", "25\u00B0S", "5\u00B0S",
                                "5\u00B0N", "25\u00B0N", "45\u00B0N", "65\u00B0N"), 
                     expand = c(0, 0),
                     name = "",
                     sec.axis = dup_axis())
setwd("~/Desktop/Postdoc/changepoint detection/modifacation")
save_plot("Figure5a_new.png", figure5a, base_height = 3, base_width = 7)

#(b)

world <- map_data("world")
max_lat = 65; max_long = 180
setwd("~/Desktop/Postdoc/changepoint detection/Main - Functional Changepoint Detection Code-selected")
load("grids_65_10_120.RData")
res = grids[select_grids, c('id','lat_lower','lat_upper','long_lower','long_upper')]
setwd("~/Desktop/Postdoc/changepoint detection/geosphere")
prec_temp = read.csv("qnormc.csv", header = TRUE)[,-1]
r_thin = seq(15001, 20000, 10)
c.mean = apply(pnorm(as.matrix(prec_temp[r_thin,])), 2, mean)
res$Changepoint = round(as.vector(c.mean)*43)
range(res$Changepoint)
p_temp = read.csv("p.csv", header = TRUE)[,-1]
p.mean = apply(as.matrix(p_temp[r_thin,]), 2, mean)
res$Type = ifelse(as.vector(p.mean) > 0.4, "alt", "null")
res$Type = ifelse(as.vector(p.mean) > 0.5, "alt", "null")
res$Type = factor(res$Type, levels = c("alt", "null"))
res$altChangepoint = res$Changepoint
tau = read.csv("tau.csv", header = TRUE)[,-1]
tau_expect <- apply(as.matrix(tau[r_thin,]), 2, mean)

res$scale = tau_expect
colnames(res)[9] <- 'temperature'
res$temperature[is.na(res$altChangepoint)] <- NA
br = seq(0,8,1)
lb = c("0","1","2","3","4","5","6","7","8")
figure5b <- ggplot() +
  coord_fixed(ratio = 1.4)+
  geom_rect(data = res, aes(xmin = long_lower, xmax = long_upper, ymin = lat_lower, ymax = lat_upper, fill = temperature))+
  #scale_fill_gradientn(colors = c("yellow", "orange", "red"), 
  #                     name = "Temperature \n Change Magnitude", na.value = "grey50", guide = "colourbar") + 
  scale_fill_gradient2(breaks = br, limits = c(0, 8), labels = lb, high = "red", mid = "yellow",midpoint = 0,
                       na.value = "grey50", name = "Temperature \n Change Magnitude") + 
  geom_map(data = world, map = world, aes(long, lat, map_id = region),
           color = "black", fill = NA, size = 0.1) +
  ylim(c(-max_lat, max_lat))+xlim(c(-max_long, max_long))+
  geom_point(aes(y=15.1383,x=120.3500),colour="red", shape=17, size = 2)+
  xlab("Longitude")+ylab("Latitude")+
  geom_text(data = grid_center[select_grids,], aes(long_lower+10, (lat_lower+lat_upper)/2, label = id, fontface=2), hjust=0.6, vjust=0.6, 
            size=2.5, col=1) + 
  theme(legend.text=element_text(size=6), 
        legend.key.height= unit(1.1, 'cm'),
        legend.title = element_text(size = 7),
        axis.title =element_text(size=8),
        axis.text.x.top = element_blank(),    # Remove top longitude labels
        axis.ticks.x.top = element_blank())+    # Remove top longitude ticks)+ 
  scale_x_continuous(limits = c(-180, 180), breaks = c(-180, -60, 60, 180), 
                     labels = c("180\u00B0W", "60\u00B0W", "60\u00B0E", "180\u00B0E"), 
                     expand = c(0, 0), 
                     name = "",
                     sec.axis = dup_axis())+
  scale_y_continuous(breaks = c(seq(-65, -5, 20),seq(5, 65, 20)), 
                     labels = c("65\u00B0S", "45\u00B0S", "25\u00B0S", "5\u00B0S",
                                "5\u00B0N", "25\u00B0N", "45\u00B0N", "65\u00B0N"), 
                     expand = c(0, 0),
                     name = "",
                     sec.axis = dup_axis())
setwd("~/Desktop/Postdoc/changepoint detection/modifacation")
save_plot("Figure5b_new.png", figure5b, base_height = 3, base_width = 7)
setwd("~/Desktop/Postdoc/changepoint detection/modifacation")
figure5 <- plot_grid(figure5a, figure5b, nrow = 2, 
                     labels = c("(a)","(b)"), label_size = 16)
save_plot("Figure5_new.png", figure5, base_height = 6, base_width = 7,bg = "white")

p1 <- ggdraw(figure5a) + 
  draw_label("(a)", x = 0.14, y = 0.87, hjust = 0, vjust = 1, size = 12)

p2 <- ggdraw(figure5b) + 
  draw_label("(b)", x = 0.10, y = 0.91, hjust = 0, vjust = 1, size = 12)

# Combine with plot_grid
figure5 <- plot_grid(p1, 
          ggdraw() + draw_plot(p2, x = 0.04, y = 0.08, width = 0.95, height = 1),  # Shift figure5b to the right
          nrow = 2,
          rel_heights = c(1, 1))  # Adjust heights if needed
save_plot("Figure5_new.png", figure5, base_height = 6, base_width = 7,bg = "white")

##############



bat_home <- batting_2024_home_pf[,c(1,5,6,9,3)]
colnames(bat_home) <- c('Team', 'home_RS', 'home_HS', 'home_HRS', 'homeG')
bat_home$Team <- gsub("[0-9]", "", bat_home$Team)
team_name_ID <- data.frame(Team = bat_home$Team, 
                           ID = c('LAD', 'NYY', 'BAL', 'PHI', 'SDP', 
                                  'HOU', 'NYM', 'MIN', 'CIN', 'CLE', 
                                  'COL', 'ARI', 'LAA', 'MIL', 'TEX', 
                                  'ATL', 'BOS', 'SEA', 'OAK', 'MIA', 
                                  'DET', 'TOR', 'PIT', 'STL', 'CHC', 
                                  'WSN', 'KCR', 'TBR', 'SFG', 'CHW'))
bat_home_ID <- merge(bat_home, team_name_ID) %>% select(-Team)
colnames(bat_home_ID)[5] <- 'teamID'

bat_away <- batting_2024_away_pf[,c(1,5,6,9,3)]
colnames(bat_away) <- c('Team', 'away_RS', 'away_HS', 'away_HRS', 'awayG')
bat_away$Team <- gsub("[0-9]", "", bat_away$Team)
bat_away_ID <- merge(bat_away, team_name_ID) %>% select(-Team)
colnames(bat_away_ID)[5] <- 'teamID'

pit_home <- pitching_2024_home_pf[,c(1,13,14,16)]
colnames(pit_home) <- c('Team', 'home_RA', 'home_HA', 'home_HRA')
pit_home$Team <- gsub("[0-9]", "", pit_home$Team)
pit_home_ID <- merge(pit_home, team_name_ID) %>% select(-Team)
colnames(pit_home_ID)[4] <- 'teamID'

pit_away <- pitching_2024_away_pf[,c(1,13,14,16)]
colnames(pit_away) <- c('Team', 'away_RA', 'away_HA', 'away_HRA')
pit_away$Team <- gsub("[0-9]", "", pit_away$Team)
pit_away_ID <- merge(pit_away, team_name_ID) %>% select(-Team)
colnames(pit_away_ID)[4] <- 'teamID'

bat_pit <- merge(bat_home_ID, merge(bat_away_ID, merge(pit_home_ID, pit_away_ID)))

bat_pit_pf <- bat_pit %>% 
  mutate(park_index_R = ((home_RS+home_RA)/(homeG)) / ((away_RS+away_RA)/(awayG)), 
         park_index_H = ((home_HS+home_HA)/(homeG)) / ((away_HS+away_HA)/(awayG)), 
         park_index_HR = ((home_HRS+home_HRA)/(homeG)) / ((away_HRS+away_HRA)/(awayG)))
n = 30
bat_pit_pf <- bat_pit_pf %>% 
  mutate(park_factor_R = ((park_index_R +1)/2)/((park_index_R+n-1)/n), 
         park_factor_H = ((park_index_H +1)/2)/((park_index_H+n-1)/n), 
         park_factor_HR = ((park_index_HR +1)/2)/((park_index_HR+n-1)/n))

park_factor_2024 <- bat_pit_pf %>% mutate(yearID = 2024, teamIDBR = teamID) %>% 
  select(yearID, teamID, park_factor_R, park_factor_H, park_factor_HR, teamIDBR)





