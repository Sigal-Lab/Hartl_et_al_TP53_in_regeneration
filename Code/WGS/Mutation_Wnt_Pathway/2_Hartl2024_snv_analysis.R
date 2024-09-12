# script to plot mouse SNVs for Hartl et al. (2024) 

library(tidyverse)
library(reshape2)
library(ggarange)
library(vcfR)

# plotting function
plot_snv = function(vcf, gene_lists, gene_drop = F, add_count = F){
  print(length(gene_lists))
  plot_input_combined = c()
  for(gene_list_name in names(gene_lists)){
    print(gene_list_name)
    gene_list = rev(sort(unlist(gene_lists[gene_list_name],use.names = F)))
    
    # 1. count exonic event per gene per sample
    plot_input = vcf %>% 
      filter(Gene.refGene %in% gene_list) %>% 
      separate(AD, into = c('RO',"AO"),sep=",") %>% 
      filter(as.numeric(AO) >= 3) %>% 
      group_by(Sample, Gene.refGene) %>% 
      mutate(ExonicFunc.refGene = ifelse(ExonicFunc.refGene=='.',Func.refGene,ExonicFunc.refGene)) %>% 
      filter(Gene.refGene %in% gene_list) %>% 
      group_by(Sample,Gene.refGene) %>%  
      dplyr::add_count(ExonicFunc.refGene) %>% 
      reshape2::dcast(Sample + Gene.refGene ~ ExonicFunc.refGene)
    
    # 2. add columns for plotting
    plot_input$plot_col = NA
    plot_input$plot_num = NA
    
    for(i in 1:nrow(plot_input)){
      mywords = plot_input[i,3:ncol(plot_input)] %>% 
        reshape2::melt() %>% 
        filter(value >= 1) %>% 
        pull(variable) %>% 
        as.character() %>% 
        paste(collapse = ',')
      plot_input$plot_col[i] = mywords
      
      mynums = plot_input[i,3:ncol(plot_input)] %>% 
        reshape2::melt() %>% 
        filter(value >= 1) %>% 
        pull(value) %>% 
        paste(collapse = ',')
      plot_input$plot_num[i] = mynums
    }
    plot_input$facet_name = gene_list_name
    plot_input_combined = rbind(plot_input_combined, plot_input)
  }
  
  # 3. plot
  gene_order = rev(unique(sort(unlist(gene_lists,use.names = F))))
  mutation_types =  c("intronic", "synonymous SNV", "nonsynonymous SNV", "nonframsehfit indel", "frameshift indel", 
                      "stop gain", "stop loss", "splicing")
  plot_input_combined$plot_col = factor(plot_input_combined$plot_col, levels = mutation_types)
  
  myplot = ggplot(plot_input_combined,
                  aes(x =factor(Sample,levels =  sort(unique(vcf$Sample))),
                      y = factor(Gene.refGene, levels = gene_order), fill = plot_col)) +
    scale_x_discrete(drop = FALSE) + 
    scale_y_discrete(drop = FALSE) + 
    scale_fill_manual(
      values = c("#6C6EA0", "#D8D174", "firebrick3", "#F0B67F", "pink", "#FFD2FC", "#9DBA8D", "#7CDEDC"),
      labels = mutation_types, 
      drop = FALSE
    ) + 
    facet_wrap(.~facet_name, scales = 'free') + 
    geom_tile(na.rm = F) +
    theme_bw() + 
    theme(legend.position = 'bottom') +
    theme(legend.spacing.y = unit(.01, 'npc'))  +
    ## important additional element
    guides(fill = guide_legend(byrow = TRUE)) + 
    labs(y = 'Gene', x = 'Sample', fill = 'Mutation')
  
  if(!gene_drop) myplot = myplot + scale_y_discrete(drop = FALSE)
  if(add_count) myplot = myplot +   geom_text(aes(label = plot_num), vjust = .5, check_overlap = T) 
  
  print(myplot)
}

# load vcf function
vcf_load = function(datafile){
  vcf_load = read.vcfR(datafile,verbose = F)
  message("vcf file loaded")
  
  if(vcfR::nrow(vcf_load)==0){
    warning("no variants detected, nothing to return")
    return(NULL)
  }else{
    infocols = sub("=.*", "",  unlist(strsplit(vcf_load@fix[1,"INFO"],";"), use.names = F))
    for(i in 0:length(infocols)){
      if(any(duplicated(infocols))){
        mydup = head(infocols[duplicated(infocols)],1)
        if(length(infocols[infocols==mydup])==2){
          infocols[tail(which(infocols==mydup),1)] = paste0(mydup,"_2")
        }else{
          stop("more than one duplicate column name found")
        }
      }else{
        message(paste("after",i,"iterations, no (more) duplicated column names found"))
        break
      }
    }
    rmnames = paste0(infocols,rep("=",length(infocols)), collapse = "|")

    fix_df = as.data.frame(getFIX(vcf_load,getINFO = T)) %>% 
      mutate(INFO = gsub(rmnames, "", INFO)) %>% 
      separate(INFO,sep = ";", into = infocols)
    
    message("fixed df made")
    
    # test if there are more than 1 samples in the .vcf file
    my_samps = colnames(vcf_load@gt)[colnames(vcf_load@gt)!="FORMAT"]
    if(length(my_samps)>1){
      multi_samps = T
      message(length(colnames(vcf_load@gt)[colnames(vcf_load@gt)!='FORMAT']), " samples detected in vcf")
    }else{
      multi_samps = F
    }
    
    mycols = unlist(strsplit(head(vcf_load@gt)[1,1],":"))
    vcf_final = list()
    
    for(mysamp in my_samps){
      message("working on ", mysamp)
      
      format_df = vcf_load@gt %>% 
        as.data.frame() %>% 
        dplyr::select(format_column = !!as.name(mysamp)) %>%
        separate(format_column, into = mycols, sep = ":")
      
      # remove one copy of a column that appears in both dataframes
      for(appears_twice in names(format_df)[names(format_df) %in% names(fix_df)]){
        print(appears_twice)
        if(!all(format_df[,appears_twice] == fix_df[,appears_twice],na.rm = T)){
          if(appears_twice %in% c("DP","AD")){
            fix_df = dplyr::select(fix_df,-c(!!as.name(appears_twice))) 
          } else {
            if(all(colnames(format_df) %in% c("GT","AD","AF","DP","F1R2","F2R1","FAD","SB"))){
              if(appears_twice == "AF") format_df = dplyr::select(format_df,-c(!!as.name(appears_twice))) 
              
            }else{
              stop(paste("ERROR: the column",appears_twice,"has differing values between the fixed and format parts of the vcf"))
            }
          }
        }else{
          message(paste("the column",appears_twice,"appears twice in the VCF, removing one instance now"))
          format_df = format_df %>% 
            dplyr::select(-c(!!as.name(appears_twice)))
        }
      }
      
      message("format df made")
      
      mylist = list("fix" = fix_df, "format" = format_df)
      vcf_combined = do.call(cbind, mylist)
      names(vcf_combined) = gsub("fix.|format.","",names(vcf_combined))
      
      if(length(my_samps)>1){
        vcf_final[[mysamp]] = vcf_combined
      }else{
        vcf_final = vcf_combined
      }
    }
    return(vcf_final)
    message("vcf load completed")
  }
}


######
## 1. define gene lists
#####

# manually curated Wnt signalling genes/genes of interest

# Oncoprint selected genes Wnt signaling
oncoprint_genes = toupper(c('Amer1','Apc','Arid1a','Axin2','Bcl9','Bcl9l','Ctnnb1','Dkk1','Dkk2','Dkk3','Dkk4','Fbxw7',
                            'Fzd10','Mir34a','Mir34b','Mir34c','Rnf43','Sox9','Tcf7l2','Trp53'))

# colitis associated CRC KIH 
colitis_genes = c('FUT8','MAGEC3','PODN','GYS1','COL6A2','BOX','CRB1','FSTL5','PCDH9','APC','LSAMP',
                  'CCSER1','FHIT','IMMP2L','MACROD2','PIBF1','OSMR','LIFR','CDKN2B','SMARCA2','SMAD4',
                  'FOXA1','NFKBIZ','PIGR','ARID1A','IL17RA','MYC','IDH1')

yap_other_genes = toupper(c('ARID1A','LKB1','FAT1','Wnt1','APC','PCP'))

yap_genes = toupper(c('YAP1','TAZ','MST1','MST2','LATS1','LATS2','NF2','Sav1','Mob1','TEAD1','RhoA',
                      'TEAD4','Ptpn14','WWTR1','RASSF1','RASSF2','RASSF5','Nore1',yap_other_genes))

######
## 2. load data
######

# load vcfs that were annotated by annovar_script_share.sh
vcf_KOWIK1_load = vcf_load('P3397_DNA_03_5001KOWIK1_pass_annotated.vcf.gz')

vcf_KOWIK1 = vcf_KOWIK1_load %>% 
    mutate(Gene.refGene = toupper(Gene.refGene),
           Sample = 'KOWIK1_mutect')


vcf_KOWIK3_load = vcf_load('P3397_DNA_04_5001KOWIK3_pass_annotated.vcf.gz')

vcf_KOWIK3 = vcf_KOWIK3_load %>% 
    mutate(Gene.refGene = toupper(Gene.refGene),
           Sample = 'KOWIK3_mutect')


######
## 3. plot
######

vcf_KOWI_combined = rbind(vcf_KOWIK3, vcf_KOWIK1)

# make plots 
plot1 = plot_snv(vcf = vcf_KOWI_combined, 
         gene_lists = list('Oncoprint genes' = oncoprint_genes), 
         add_count = F)

plot2 = plot_snv(vcf = vcf_KOWI_combined, 
                 gene_lists = list('YAP genes' = yap_genes), 
                 add_count = F)

ggpubr::ggarrange(plot1, plot2, common.legend = T, legend = 'bottom')

# save plots
png(filename = 'mouse_snv_plot.png',
    width = 2000, height = 1500, res = 250)
ggpubr::ggarrange(plot1, plot2, common.legend = T, legend = 'bottom')
dev.off()


plot_snv(vcf = vcf_KOWI_combined, gene_list = oncoprint_genes, add_count = T)
plot_snv(vcf = vcf_KOWI_combined, gene_list = list('colitis_genes'=colitis_genes), add_count = T)
plot_snv(vcf = vcf_KOWI_combined, gene_list = list('yap' = yap_genes), add_count = T)


