library(hexSticker)
library(sismonr)
library(tidyverse)
library(scales)
library(patchwork)
library(wesanderson)
library(ggnetwork)


#############################################
##          SIMULATE DATA TO PLOT          ##
#############################################

set.seed(152)

grn <- createInSilicoSystem(empty = TRUE,
                            G = 3,
                            PC.p = 1,
                            PC.TC.p = 1)

grn$genes$TCrate <- c(1, 0.005, 0.005)
grn$genes$RDrate <- c(0.02, 0.01, 0.01)
grn$genes$TLrate <- c(0.1, 0.1, 0.1)
grn$genes$PDrate <- c(0.01, 0.01, 0.01)

grn <- addEdge(grn, "1", "2", regsign = "1", kinetics = list(TCfoldchange = 100,
                                                             TCbindingrate = 0.1,
                                                             TCunbindingrate = 0.05))
grn <- addEdge(grn, "1", "3", regsign = "1", kinetics = list(TCfoldchange = 30,
                                                             TCbindingrate = 0.1,
                                                             TCunbindingrate = 0.05))
grn <- addEdge(grn, "2", "3", regsign = "1", kinetics = list(TCfoldchange = 30,
                                                             TCbindingrate = 0.1,
                                                             TCunbindingrate = 0.05))
grn <- addEdge(grn, "3", "1", regsign = "-1", kinetics = list(TCfoldchange = 0,
                                                              TCbindingrate = 0.5,
                                                              TCunbindingrate = 0.05))


pop <- createInSilicoPopulation(1,
                                grn,
                                InitVar = list(R = rep(0, 3),
                                               P = rep(0, 3)))

sim <- simulateInSilicoSystem(grn,
                              pop,
                              50,
                              ntrials = 100)

#############################################
##             PLOT NETWORK                ##
#############################################

border_cols <- "black"

set.seed(2)
p_grn <- ggnetwork(getEdges(grn), arrow.gap = 0.12) %>%
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(linetype = RegSign),
             curvature = 0.2,
             arrow = arrow(type = "closed", length = unit(3, "pt")),
             size = 0.2,
             colour = border_cols) +
  geom_nodes(aes(fill = vertex.names),
             size = 2.5,
             shape = 21,
             colour = border_cols,
             stroke = 0.2) +
  scale_fill_viridis_d(option = "plasma", direction = -1, guide = "none") +
  scale_linetype_manual(values = c("1" = 1, "-1" = 2), guide = "none") +
  scale_x_continuous(limits = c(-0.5, 1.5)) +
  scale_y_continuous(limits = c(-0.4, 1.4)) +
  theme_void() +
  theme_transparent()

#############################################
##           PLOT SIMULATION               ##
#############################################

sim_df <- mergeAlleleAbundance(sim$Simulation) %>%
  as_tibble() %>%
  pivot_longer(cols = matches("^(R|P)"),
               names_to = "mol",
               values_to = "abundance") %>%
  mutate(abundance = abundance + 0.1,
         gene = str_extract(mol, "\\d"),
         type = case_when(str_detect(mol, "R") ~ "mRNAs",
                          TRUE ~ "Proteins")) %>%
  group_by(time, gene, type) %>%
  summarise(mean = mean(abundance),
            lq = quantile(abundance, probs = 0.1),
            hq = quantile(abundance, probs = 0.9))

make_plot <- function(what, topmargin){
  sim_df %>%
    filter(type == what) %>%
    ggplot(aes(x = time)) +
    geom_ribbon(aes(ymin = lq, ymax = hq, fill = gene), alpha = 0.4) +
    geom_line(aes(y = mean, colour = gene)) +
    facet_grid(type ~ ., scales = "free_y") +
    #scale_colour_brewer(palette = "Set1", guide = "none") +
    #scale_fill_brewer(palette = "Set1", guide = "none") +
    scale_colour_viridis_d(option = "plasma", direction = -1, guide = "none") +
    scale_fill_viridis_d(option = "plasma", direction = -1, guide = "none") +
    #scale_colour_manual(values = wes_palette("Royal1", 3, "discrete"), guide = "none") +
    #scale_fill_manual(values = wes_palette("Royal1", 3, "discrete"), guide = "none") +
    scale_x_continuous(expand = expansion()) +
    scale_y_log10(labels = label_comma(drop0trailing = TRUE),
                  expand = expansion()) +
    theme_void() +
    theme(axis.line = element_line(colour = border_cols, size = 0.2, arrow = arrow(type = "closed", length = unit(2, "pt"))),
          strip.background = element_blank(),
          strip.text = element_blank(),
          plot.margin = unit(c(topmargin, 0, 0, 0), "pt"))
}

p <- p_grn +
  (make_plot("Proteins", 0) /
     make_plot("mRNAs", 5))


name_col <- "#002966"
sticker(
  p & theme_transparent(),
  package = "sismonr",
  p_size = 20,
  p_color = name_col,
  s_x = 0.9,
  s_y = 0.85,
  s_width = 1.8,
  s_height = 0.8,
  h_color = name_col, # "#e6f0ff",
  h_fill = "#cce0ff", # "#66a3ff",
  filename = "inst/figures/sismonr.png"
)




