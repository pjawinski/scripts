#!/usr/bin/env Rscript
# ======================================================
# === plot height against polygenic score for height ===
# ======================================================
# feel free to adjust this script for your own purposes

message('--- plot height vs. polygenic score ---')

# load required packages
for (pkg in c('dplyr','ggplot2','gganimate')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# load data
height = read.delim('data/phenotypes/height.txt', header = T) # text file with columns IID and height
pgs = read.delim('results/height/pgs/pgs.score', header = T) # output file created with pgsheight.sh
sex = read.delim('data/genetics/processed/imp_qc/chr1.psam', header = T) # file containing columns IID and sex

# remove FID column
pgs = pgs[,c('IID','CHR_COUNT','NMISS_ALLELE_CT','SCORE1_SUM')]
sex = sex[,c('IID','SEX')]

# merge datasets
df = height %>% left_join(sex, by = 'IID') %>% left_join(pgs, by = 'IID')
df = df[complete.cases(df),]

# residualise by sex
df$heightResiduals = mean(df$Phenotype) + residuals(lm(df$Phenotype ~ df$SEX))
df$pgsResiduals = scale(residuals(lm(df$SCORE1_SUM ~ df$SEX)))

# calculate correlation
R2 = cor(df$heightResiduals,df$pgsResiduals)^2

# draw plot
message('Plotting pgs vs. phenotype (png).')
pl = ggplot(df, aes(x = pgsResiduals, y = heightResiduals)) +
  geom_point(colour = '#4393C3', alpha = 1, size = 1) + # #13848F
  geom_smooth(method = mgcv::gam, formula = y ~ x, se = TRUE, linewidth = 0.25, color = "black", alpha = 0.2) +
  scale_x_continuous(limits = c(-3.5,3), breaks = seq(-3,3,2)) +
  scale_y_continuous(limits = c(155,195), breaks = seq(150,190,10)) +
  theme_bw(base_size=10) +
  theme(plot.margin = unit(c(5, 5, 5, 5), "mm"),
        panel.border = element_rect(linewidth = 0.25),
        line = element_line(linewidth = 0.25),
        plot.title = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 1, lineheight = 1.1, face ='bold'),
        axis.title.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = -0.5, lineheight = 1.1),
        axis.title.y = element_text(size = 12, angle = 90, hjust = 0.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks = element_line(linewidth = 0.25)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.5, vjust = 4, label = paste0('R2 = ',sprintf('%.3f', R2)), size = 3) +
    labs(title = 'Genetic prediction of height', 
       x = 'Polygenic score (z-scale)',
       y = 'Height (cm)')

# save file
message('Writing .png file.')
png(width = 4.0, height = 3.2, units = "in", res = 300, filename = 'results/height/pgs/pgs.png'); pl; invisible(dev.off())

# create animation
message('Preparing data for animation.')
for (i in c(1:22)) {
  tmp = read.delim(paste0('results/height/pgs/chr/chr',i,'.sscore'), header = T)
  tmp = tmp[,c('IID','NMISS_ALLELE_CT','SCORE1_SUM')]
  tmp = height %>% left_join(sex, by = 'IID') %>% left_join(tmp, by = 'IID')
  tmp = tmp[complete.cases(tmp),]
  tmp = tmp[,-which(names(tmp) %in% 'FID')]
  tmp$CHR_COUNT = i
  
  if (i == 1) {
    df = tmp 
    last = tmp
  } else { 
    for (id in tmp$IID) {
      tmp$SCORE1_SUM[tmp$IID == id] = last$SCORE1_SUM[last$IID == id] + tmp$SCORE1_SUM[tmp$IID == id]
    }
    df = rbind(df,tmp)
    last = tmp
  }
}
  
# z-transform values
df$R2 = ""
for (i in unique(df$CHR_COUNT)) {
  df$heightResiduals[df$CHR_COUNT == i] = mean(df$Phenotype[df$CHR_COUNT == i]) + residuals(lm(df$Phenotype[df$CHR_COUNT == i] ~ df$SEX[df$CHR_COUNT == i]))
  df$pgsResiduals[df$CHR_COUNT == i] = residuals(lm(df$SCORE1_SUM[df$CHR_COUNT == i] ~ df$SEX[df$CHR_COUNT == i]))
  df$pgsResiduals[df$CHR_COUNT == i] = scale(df$pgsResiduals[df$CHR_COUNT == i])
  df$R2[df$CHR_COUNT == i & !duplicated(df$CHR_COUNT)] = paste0("R2 = ", sprintf('%.3f',cor(df$heightResiduals[df$CHR_COUNT == i],df$pgsResiduals[df$CHR_COUNT == i])^2))
}

# creating animated plot
message('Plotting pgs vs. phenotype (gif).')
pl = ggplot(df, aes(x = pgsResiduals, y = heightResiduals, label = R2)) +
  geom_point(colour = '#4393C3', alpha = 0.8, size = 1) +
  scale_x_continuous(limits = c(-3.5,3), breaks = seq(-3,3,2)) +
  scale_y_continuous(limits = c(155,195), breaks = seq(150,190,10)) +
  theme_bw(base_size=10) +
  theme(plot.margin = unit(c(5, 5, 5, 5), "mm"),
        panel.border = element_rect(linewidth = 0.25),
        line = element_line(linewidth = 0.25),
        plot.title = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 1, lineheight = 1.1, face ='bold'),
        axis.title.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = -0.5, lineheight = 1.1),
        axis.title.y = element_text(size = 12, angle = 90, hjust = 0.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks = element_line(linewidth = 0.25)) +
  geom_text(aes(-3, 190), size = 3.5, hjust = 0, color = "black") +
  transition_time(CHR_COUNT) +
  ease_aes('cubic-in-out') +
  labs(title = 'Genetic prediction of height', 
       x = 'Polygenic score (chr1-{frame_time})',
       y = 'Height (cm)')

# animate(pl, width = 4.0, height = 3.5, units = "in", res = 60, end_pause = 30, duration = 5, fps = 10)
message('Creating animation and writing .gif file.')
pl = animate(pl, width = 4.0, height = 3.5, units = "in", res = 200, end_pause = 30, duration = 20, fps = 10)
anim_save('results/height/pgs/pgs.gif', pl)
message('--- Completed: plot height vs. polygenic score ---')

  


