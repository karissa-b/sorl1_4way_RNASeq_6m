library(ggpubr)
library(tidyverse)
library(magrittr)
library(ggrepel)
library(UpSetR)

chrst <- topTables_RUVk1 %>% 
  lapply(function(x) {
    x %>% 
      left_join(grGenes %>% as_tibble() %>% dplyr::select(gene_id, chromosome = seqnames, start, end)
      ) %>% 
      dplyr::filter(chromosome %in% seq(1:25)) %>% 
      group_by(gene_id) %>% 
      mutate(mid = mean(c(start, end))) %>% 
      ungroup %>% 
      mutate(chromosome = factor(chromosome, levels = seq(1:25)))  %>% 
      group_by(chromosome) %>% 
      summarise(chrLen = max(mid)) %>% 
      mutate(chrSt = cumsum(chrLen)-chrLen) %>% 
      dplyr::select(-chrLen) 
  }) %>% 
  bind_rows() %>% 
  unique

man <- topTables_RUVk1 %>% 
  lapply(function(x) {
    x %>% 
      left_join(grGenes %>% as_tibble() %>% dplyr::select(gene_id, chromosome = seqnames, start, end)
      ) %>% 
      dplyr::filter(chromosome %in% seq(1:25)) %>% 
      left_join(chrst %>% 
                  mutate(chromosome = factor(chromosome, levels = seq(1:25))) 
      ) %>% 
      group_by(gene_id) %>% 
      mutate(mid = mean(c(start, end))) %>% 
      ungroup %>% 
      dplyr::arrange(mid, chromosome) %>% 
      mutate(midCum = chrSt + mid)
  })


axis <- 
  man$`R122Pfs/+` %>%
  mutate(chromosome = factor(chromosome, levels = seq(1:25))) %>% 
  group_by(chromosome) %>%
  dplyr::arrange(chromosome) %>%
  summarize(center = (max(midCum) + min(midCum)) / 2) %>%
  mutate(
    colour = rep(c("grey40", "black"), length.out = 25),
    colour = ifelse(chromosome == 15, "red", colour)
  )

man %>% 
  bind_rows() %>%  
  mutate(chromosome = factor(chromosome, levels = seq(1:25))) %>% 
  ggplot(aes(x = midCum, y = -log10(PValue))) +
  theme_bw() +
  geom_point(aes(color=chromosome), alpha = 0.5, size = 1) +
  scale_color_manual(values = axis$colour ) +
  scale_x_continuous(label = axis$chromosome, breaks = axis$center) +
  coord_cartesian(ylim = c(0,15)) +
  labs(x = "Chromosome", y = expression(paste(-log[10], "(p)"))) +
  facet_wrap(~coef) +
  geom_text_repel(data = dplyr::filter(bind_rows(man), DE == TRUE),
                  aes(label = symbol)
  ) +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  ggsave("manhatplot_sorl1.png", width = 30, height = 10, units = "cm", dpi = 900, scale = 1.5)



ploteofad <- ggplot(man_eofad, aes(x = midCum, y = -log10(PValue))) +
  theme_bw()+
  geom_point(aes(color=chromosome), alpha = 0.5, size = 1) +
  scale_color_manual(values = axis_eofad$colour) + 
  scale_x_continuous(label = axis_eofad$chromosome, breaks = axis_eofad$center) +
  coord_cartesian(ylim = c(0,15)) +
  labs(x = "Chromosome", y = expression(paste(-log[10], "(p)"))) +
  geom_hline(yintercept = -log10(2.82e-5)) +
  geom_text_repel(data = dplyr::filter(man_eofad, DE), 
                  aes(label = gene_name)
  ) +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  ggtitle("EOfAD-like") 

ggarrange(plotfai, ploteofad, ncol = 1) +
  ggsave("stackedmanhats.png", width = 20, height = 20, units = "cm", dpi = 900, scale = 2)
