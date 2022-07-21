# # Plot libraries
# library(ggrepel)
# library(broom)
# library(gridExtra)
# library(ggplot2)
# library(scales)
# library(grid)
# library(colorspace)

# library(data.table)
# library(stringr)
# library(dplyr)
# library(survival)
# library(survminer)
# library(logistf)
# library(base)
# library(cmprsk)
# library(prodlim)
# library(bshazard)
# library(mstate)

# library(doMC)
# library(glmnet)

# # Sankey Plot libraries
# library(rbokeh)
# library(networkD3)

### Main helper file ###


df_final <- read.table("data/paper_aml_prognosis_updated.tsv",check.names = FALSE)
df_final_bis <- read.table("data/aml_prognosis_updated.tsv")

df_validation <- read.table("data/paper_full_data_validation.tsv",check.names = FALSE)
df_validation_bis <- read.table("data/full_data_validation.tsv")




genes <- colnames(df_final)[grep("^[A-Z]", colnames(df_final))]
cytos <- c(colnames(df_final)
           [c(grep("^t\\(", colnames(df_final)),
            grep("^\\+", colnames(df_final)),
            grep("^\\-", colnames(df_final)),
            grep("^inv\\(+", colnames(df_final)))],
            "complex","others_transloc")

cytos_bis <- c(colnames(df_final_bis)
           [c(grep("^t_", colnames(df_final_bis)),
            grep("^add_", colnames(df_final_bis)),
            grep("^del_", colnames(df_final_bis)),
            grep("^minus", colnames(df_final_bis)),  
            grep("^inv_", colnames(df_final_bis)))],
            "complex","others_transloc")

eln <- colnames(df_final)[grep("^eln_2017_", colnames(df_final))]
# comp <- colnames(df_final)[c(grep("^full_component_", colnames(df_final)),"overlap")]

princ_comps <- colnames(df_final)[grep("^principal_component_", colnames(df_final))]
princ_comps_bis <- colnames(df_final_bis)[grep("^principal_component_", colnames(df_final_bis))]

clin <- c("ahd","perf_status","bm_blasts","secondary","wbc","hb","plt")      
demo <- c("gender","age")

vect <- apply(X=df_final[,genes],2,FUN=function(x) 100*length(which(x==1))/dim(df_final)[1])
gen <- colnames(df_final[,match(names(vect[vect>=2]),names(df_final))])
              
vect <- apply(X=df_final[,cytos],2,FUN=function(x) 100*length(which(x==1))/dim(df_final)[1])
cyto <- colnames(df_final[,match(names(vect[vect>=2]),names(df_final))])
              
vect <- apply(X=df_final_bis[,cytos_bis],2,FUN=function(x) 100*length(which(x==1))/dim(df_final_bis)[1])
cyto_bis <- colnames(df_final_bis[,match(names(vect[vect>=2]),names(df_final_bis))])




cytos_val <- c(colnames(df_validation)
           [c(grep("^t\\(", colnames(df_validation)),
            grep("^\\+", colnames(df_validation)),
            grep("^\\-", colnames(df_validation)),
            grep("^inv\\(+", colnames(df_validation)))],
            "complex")
              
cytos_val_bis <- c(colnames(df_validation_bis)
           [c(grep("^t_", colnames(df_validation_bis)),
            grep("^add_", colnames(df_validation_bis)),
            grep("^del_", colnames(df_validation_bis)),
            grep("^minus", colnames(df_validation_bis)),  
            grep("^inv_", colnames(df_validation_bis)))],
            "complex")

# cytos_val
genes_val <- colnames(df_validation)[grep("^[A-Z]", colnames(df_validation))]
clin_val <- c("BM_Blasts","HB","perf_status","PLT","WBC","secondary","ahd")      
demo_val <- c("Gender","Age")
genes_val <- genes_val[!genes_val %in% c(clin_val,demo_val,"OS","OS_Status")] # remove uppercases from clin,demo,OS,OS_Status

color_values <- list()
color_values["Trisomies"] <- "olivedrab3"
color_values["CEBPA_bi"] <- "darkorange"
color_values["sAML1"] <- "#009E73"
color_values["sAML2"] <- "#0072B2"                  
color_values["DNMT3A_IDH1_2"] <- "#CC79A7"
color_values["inv(16)"] <- "cornflowerblue"
# color_values["inv_16"] <- "cornflowerblue"
color_values["inv(3)"] <- "mediumpurple4"
# color_values["inv_3"] <- "mediumpurple4"
color_values["no_events"] <- "tan"
color_values["mNOS"] <- "slategray"
color_values["NPM1"] <- "lavenderblush2"
# color_values["NPM1 no ITD"] <- "lavenderblush2"
# color_values["NPM1 + ITD"] <- "lavenderblush2"
color_values["t(11)"] <- "firebrick1"
# color_values["t_11"] <- "firebrick1"
color_values["t(15;17)"] <- "black"
# color_values["t_15_17"] <- "black"
color_values["t(6;9)"] <-"lightpink4"
# color_values["t(6;9) no ITD"] <-"lightpink4"
# color_values["t(6;9) + ITD"] <-"lightpink4"
# color_values["t_6_9"] <-"lightpink4" 
color_values["t(8;21)"] <- "#a6bddb"
# color_values["t_8_21"] <- "#a6bddb"
color_values["TP53_complex"] <- "firebrick4"
color_values["WT1"] <- "#fdbb84"
# color_values["WT1 no ITD"] <- "#fdbb84"
# color_values["WT1 + ITD"] <- "#fdbb84"
color_values["overlap"] <- "#870C14"
# color_values[c("_adverse","_intermediate","_favorable")] <- c("gray0","gray1","gray2")
color_values["ELN Favorable"] <- "#2b8cbe"
color_values["ELN Intermediate"] <- "#59A08B"
color_values["ELN Adverse"] <- "#EE9937"


mrd <- c("CR_MRD_neg","CR_MRD_pos","all_others")

colors_analysis <- c(eln="#EE9937",comp="#59A08B",princ_comps="#59A08B",comp_bis="#59A08B",genes="#BFBFBF",genes_bis="#BFBFBF",cytos="#2b8cbe",cytos_bis="#2b8cbe",clin="#870C14",
                     demo="#a6bddb",age="#a6bddb",gen_cyto="pink",eln_gen_cyto="#fdbb84",comp_gen_cyto="lightgoldenrod",all="lightsalmon",mrd="pink",mrdneg="#8B5F65",mrdpos="#F08080")



# Util Functions :

get_table <- function(data, add_total_count = TRUE) {
    # Print a custom sorted table (count and frequency) of a categorical feature
    # â†’ Arguments
    #     - data           : one dimensional vector
   
    #     - add_total_count: if TRUE, add a final row showing the total count

    result <- data.frame(table(data))
    colnames(result) <- c('values', 'count')
    result$freq <- result$count * 100 / sum(result$count)
    
    result <- result[order(result$count, decreasing =T),]
    if (add_total_count){
        result <- rbind(result, data.frame('values' = '-- total --',
                                            'count'  = sum(result$count),
                                            'freq'   = '100%'))
    }
    

    return (result)
}


set_notebook_plot_size <- function(width = 6, height = 3) {
    # Set the plot size in jupyter notebook
    # â†’ Arguments
    #     - width
    #     - height

    options(repr.plot.width = width, repr.plot.height = height)
}

wrap_cox <- function(columns){
    form <- paste("Surv(os, os_status) ~ ",columns[1], sep=" ")
    for (co in columns[2:length(columns)]){
        form <- paste (form,"+",co)
    }
    return(form)
}

is_one_to_one_mapping <- function(DF) {
   all(ave(DF$id, DF$patient_id, FUN = function(x) length(unique(x))) == 1)
}
           
####################################### New Proposal Classification #######################################
           
molecular_proposal_classification_without_ITD_shift <- function(df){
    
    df$molecular_classification <- "none"

    df[(df$principal_component_NPM1==1 |df$principal_component_inv_16==1 |
                df$principal_component_t_8_21==1 | df$principal_component_t_15_17==1  |
                df$principal_component_CEBPA_bi==1 | df$principal_component_no_events==1 ),"molecular_classification"] <- "favorable"

    df[(df$principal_component_sAML1==1 | df$principal_component_t_6_9==1 | df$principal_component_DNMT3A_IDH1_2==1 |
               df$principal_component_Trisomies==1 | df$principal_component_t_11==1 | df$principal_component_WT1==1 |
                df$principal_component_DNMT3A_IDH1_2==1 | df$principal_component_mNOS==1),"molecular_classification"] <- "intermediate"

    df[( df$principal_component_sAML2==1 | df$principal_component_TP53_complex==1   |
                df$principal_component_inv_3==1) ,"molecular_classification"] <- "adverse"
    
    return(df)

}
           
molecular_proposal_classification <- function(df){
    
    df$molecular_classification <- "none"

    df[(df$principal_component_NPM1==1 |df$principal_component_inv_16==1 |
                df$principal_component_t_8_21==1 | df$principal_component_t_15_17==1  |
                df$principal_component_CEBPA_bi==1 | df$principal_component_no_events==1 ),"molecular_classification"] <- "NEW_favorable"

    df[(df$principal_component_sAML1==1 | df$principal_component_t_6_9==1 | df$principal_component_DNMT3A_IDH1_2==1 |
               df$principal_component_Trisomies==1 | df$principal_component_t_11==1 | df$principal_component_WT1==1 |
                df$principal_component_DNMT3A_IDH1_2==1 | df$principal_component_mNOS==1),"molecular_classification"] <- "NEW_intermediate"

    df[( df$principal_component_sAML2==1 | df$principal_component_TP53_complex==1   |
                df$principal_component_inv_3==1) ,"molecular_classification"] <- "NEW_adverse"
    
    df[df$molecular_classification=="NEW_intermediate" & df$ITD==1,"molecular_classification"] <- "NEW_adverse"

    df[df$molecular_classification=="NEW_favorable" & df$principal_component_NPM1==1 & df$ITD==1,"molecular_classification"] <- "NEW_intermediate"

    
    return(df)

}           
           
####################################### End New Proposal Classification #######################################
           
####################################### Survival Plots + pvalues and HR #######################################
              
# Surv plot default colors
val <-c("#EE9937","#5C5C5C","#870C14","#BFBFBF","#59A08B","#2b8cbe","#a6bddb","#fdbb84","#e79f00","#000000","darkseagreen","lightskyblue","#0072B2","pink","blue","green")
              
plot_surv_curves <-function(fit,submain="",vals=val,legend="top",risk.tab=F,y="",linetype=1,size=4,pval=TRUE,pval.coord=c(0,0.05),font.legend=28,xlab="Time (years)",legend.title="",
                            risk.table.fontsize=10,font.x=24,font.y=24,font.tickslab=26,...){
    
    # Create the vector to use for leg.labs
    vec <- NULL
    for (i in 1:length(fit$strata)){
    vec <- cbind(vec,paste(str_remove(names(fit$strata[i]),"comparison="),fit$n[i],sep=" , N="))
    }
    
    ggsurvplot(fit,  pval = pval,main = "Survival curve",risk.table=risk.tab,submain = submain,palette=vals,legend=legend,pval.size=8,pval.coord=pval.coord,risk.table.fontsize=risk.table.fontsize,xlab=xlab,legend.labs=vec,ylab=y,legend.title=legend.title,...,
               linetype=linetype,size=size,
                  ggtheme = theme_survminer(
                 font.main = c(30, "plain", "black"),
                 font.submain = c(30, "plain", "black"),
                 font.legend=font.legend,
                 font.caption = c(30, "plain", "black"),
                 font.x = c(font.x, "plain", "black"),
                 font.y = c(font.y, "plain", "black"),
                 font.tickslab = c(font.tickslab, "plain", "black"),
                      base_size=30)) 

}
              
              
get_pvalue_and_HR <- function(formula,data,text="p"){
    res <- pairwise_survdiff(formula,data)
    try(res <- ifelse(round(res$p.value[1],4)==0,paste(text,"<0.0001",sep=""),paste(text,"=",round(res$p.value[1],4),sep="")))
    sum <- summary(coxph(formula,data, method = "breslow"))
    index <- ifelse(!is.na(sum$conf.int[1,1]),1,2)
    text_HR <- paste("; HR:",round(sum$conf.int[index,1],2)," (",round(sum$conf.int[index,3],2),"-",round(sum$conf.int[index,4],2),")",sep="")
    res <- paste(res, text_HR)
    return(res)
}

           
get_tdc_Transplant_pvalues <- function(formula,data){
    sum <- summary(coxph(formula = formula, data = data))
    res <- sum$coefficients["Transplant",5]
    res <- ifelse(round(res,4)==0,"p<0.0001",paste("p","=",round(res,4),""))
    text_HR <- paste("; HR :",round(sum$conf.int["Transplant",1],2)," ( ",round(sum$conf.int["Transplant",3],2),"-",round(sum$conf.int["Transplant",4],2),")")
    return(paste(res, text_HR))
}
####################################### End Survival Plots + pvalues and HR #######################################

              
####################################### Cumulative Incidence Plots #######################################              
gg_competingrisks.cuminc <- function(fit, gnames = NULL, gsep=" ",
                                    coef = 1.96, conf.int = FALSE, line.size=4, group.levels=NULL,title="Relapse risk") {
    
  # -- fit is from cuminc --  #
  if (!is.null(fit$Tests))
    fit <- fit[names(fit) != "Tests"]
  fit2 <- lapply(fit, `[`, 1:3)
  if (is.null(gnames)) gnames <- names(fit2)
  fit2_list <- lapply(seq_along(gnames), function(ind) {
    df <- as.data.frame(fit2[[ind]])
    df$name <- gnames[ind]
    df
  })
  time <- est <- event <- group <- NULL
  df <- do.call(rbind, fit2_list)
  df$event <- sapply(strsplit(df$name, split=gsep), `[`, 2)
  df$group <- sapply(strsplit(df$name, split=gsep), `[`, 1)
  df$std <- std <- sqrt(df$var)
  if (!is.na(group.levels)) {
      df$group = factor(df$group, levels=group.levels)
  }

  df <- df[df$event=="relapse",]

  pl <-  ggplot(df, aes(time, est, color=group)) + geom_line(size=line.size,linetype=1) + ylab("Cumulative incidence") + xlab("Time from CR (years)") + 
    theme_survminer(
                 font.main = c(30, "plain", "black"),
                 font.submain = c(30, "plain", "black"),
                 font.legend=28,
                 font.caption = c(30, "plain", "black"),
                 font.x = c(24, "plain", "black"),
                 font.y = c(24, "plain", "black"),
                 font.tickslab = c(26, "plain", "black"),
                      base_size=30)  +
    ggtitle(title) + theme(plot.title = element_text(hjust = 0.5,size=30)) + scale_x_continuous(breaks=seq(0, 11, 3))

  if (conf.int) {
    pl <- pl + geom_ribbon(aes(ymin = est - coef*std, ymax=est + coef*std, fill = event), alpha = 0.2, linetype=0)
  }

  return(pl)
}

####################################### End Cumulative Incidence Plots #######################################              
              
####################################### Nice Forest Plots #######################################     

pval_to_signif_code_Ref <- function (p_vals) 
{
    signifs <- ifelse(p_vals < 0.0001, "****", ifelse(p_vals < 0.001, 
        "***", ifelse(p_vals < 0.01, "**", ifelse(p_vals < 0.05, "*", "ns"))))
    
    
    return(replace(signifs,is.na(signifs),"Ref"))
}              
nice_forest_plot <- function(fit,data,colors,title){
    sum <- summary(fit)
    forest_df <- as.data.frame(cbind(sum$coef[,c(2,5)],sum$conf.int[,c(3,4)]))
    forest_df$feature <- factor(rownames(forest_df),levels=rev(unique(rownames(forest_df))))
    colnames(forest_df) <- c("HR","pval","HRlow","HRhigh","feature")
    forest_df$annot <- sapply(forest_df$pval, pval_to_signif_code_Ref)

    forest_plot <- ggplot(forest_df, aes(x=feature,y=HR,colour=feature)) + 
    geom_point(size=7) + 
    geom_hline(yintercept=1,color="darkgrey",size=1.2) + 
    geom_hline(yintercept=0.5,color="lightgrey",size=0.6) +
    geom_hline(yintercept=2,color="lightgrey",size=0.6) + 
    geom_hline(yintercept=3,color="lightgrey",size=0.6) +
    geom_errorbar(aes(ymin=HRlow,ymax=HRhigh,colour=feature),size=1.5,width=0.6) + 
    geom_text(aes(x=feature,y=0.16,label=annot),col="black",size=8) + 
    scale_y_log10(breaks=c(0.5,1,2,5,10),labels=c(0.5,1,2,5,10)) + 
    scale_colour_manual(values=colors) +
    coord_flip()+xlab("Features")+ylab(title)+
    theme_classic() + theme(text=element_text(size=25,colour="black"),
                                 axis.text=element_text(size=25,colour="black"),
                                 axis.title=element_text(size=25,colour="black"),
                                 strip.text.x = element_text(size=25),
                                 legend.position="none")
    
    return(forest_plot)
}


####################################### End Nice Forest Plots #######################################   
              
              
overall_strat_frequency <- function(df,top=70,features=c(genes,cytos)){
    tmp <- df
    df_to_plot <- tmp[,c(features,"eln_2017")] %>% group_by(eln_2017) %>% summarise_all(list(sum))
    df_to_plot$eln_2017 <- NULL
    df_to_plot <- df_to_plot[,names(sort(colSums(df_to_plot), decreasing = TRUE))]
    df_to_plot <- cbind(eln = c("Adverse","Intermediate","Favorable"),df_to_plot[,1:top])
    df_to_plot [1,2:ncol(df_to_plot)] <- 100 * df_to_plot [1,2:ncol(df_to_plot)] /nrow(tmp)
    df_to_plot [2,2:ncol(df_to_plot)] <- 100 * df_to_plot [2,2:ncol(df_to_plot)] / nrow(tmp)
    df_to_plot [3,2:ncol(df_to_plot)] <- 100 * df_to_plot [3,2:ncol(df_to_plot)] / nrow(tmp)
    gg <- melt(setDT(df_to_plot),id="eln",value.name="Frequency", variable.name="Variable")

    ggplot(gg, aes(x = Variable, y = Frequency,fill=eln)) +
        geom_bar(stat='identity')+scale_fill_manual(values=c("#EE9937","#2b8cbe","#59A08B"))+theme_bw()+theme(axis.title=element_text(hjust = 0.5,size=45,face="bold"),axis.text.x=element_text(angle=90, hjust=1),plot.title = element_text(hjust=0.5,size=55),
axis.text=element_text(size=30,face="bold"),legend.text=element_text(size=30),legend.title=element_text(size=45,face="bold"))+labs(fill="ELN 2017")+
    xlab("Genetic lesions")
    }
              
              
####################################### Univariate and multivariate regression ####################################### 
              
univariate_volcano <- function(df_final,target,features,type="continuous",quantile=c(0,1),p_value_threshold=0.05,size_title=25,legend.position="none",
                               colors_analysis=colors_analysis,change_name=T,ylab="-log10(FDR)"){
    
    df <- data.frame(beta = double(),pvalue = double(),Frequency = double()) 
    
    #     png(paste(target,"_univariate.png",sep=""),width=4500,height=2500,res=300)
    data <- df_final
    for (col in features){
        if (type=="continuous"){
            fit <- lm(as.formula(paste(paste(target," ~ ",sep= ""),col,sep="")), data=data)
            df[col,1:3] <- c(coef(fit)[[col]],glance(fit)[["p.value"]],100*sum(data[,col])/dim(data)[1])
        }else{
            fit <- logistf(formula = as.formula(paste(paste(target," ~ ",sep= ""),col,sep="")), data = data)
            df[col,1:3] <- c(fit$coef[[2]],fit$prob[[2]],100*sum(data[,col])/dim(data)[1])
            }
        }
    

    df[,"adjusted_pvalue"] <- p.adjust(df$pvalue)
    df[,"FDR"] <- -log(p.adjust(df$pvalue),10)
    df["names"] <- rownames(df)
    if(change_name){
        df$names <- ifelse(grepl("^add",df$names),str_replace(df$names,"add_","+"),
                          ifelse(grepl("^del",df$names),str_replace(df$names,"del_","-"),
                                ifelse(grepl("minusy",df$names),str_replace(df$names,"minusy","-Y"),
                                      ifelse(grepl("^t_",df$names),paste(str_replace(str_replace(df$names,"t_","t("),"_",";"),")",sep=""),
                                          ifelse(grepl("^inv_",df$names),paste(str_replace(df$names,"inv_","inv("),")",sep="")
                                                 ,df$names)))))
        
        # same for comp :
        
        df$names <- ifelse(df$names=="full_component_inv_16","full_component_inv(16)",
                           ifelse(df$names=="full_component_inv_3","full_component_inv(3)",
                                  ifelse(df$names=="full_component_t_11","full_component_t(11)",
                                         ifelse(df$names=="full_component_t_15_17","full_component_t(15;17)",
                                                ifelse(df$names=="full_component_t_6_9","full_component_t(6;9)",
                                                       ifelse(df$names=="full_component_t_8_21","full_component_t(8;21)",df$names))))))

    }

    
# change some comp names :


    df$Model <- ifelse(df$names %in% genes,"genes",
                      ifelse(df$names %in% c(cytos,cytos_bis),"cytos",
                            ifelse(df$names %in% clin, "clin",
                                  ifelse(df$names %in% demo, "demo",
                                        ifelse(df$names %in% eln, "eln",
                                               ifelse(df$names %in% mrd,"mrd","comp"))))))
    df$Frequency <- ifelse(df$names %in% clin, 5,
                                  ifelse(df$names %in% demo, 5, df$Frequency))
    df$names <- str_remove(df$names,"full_component_")
    df <- df[order(df$adjusted_pvalue),]

    s <- ggplot(df, aes(x=beta, y=FDR)) + #volcanoplot with log2Foldchange versus pvalue
            geom_point(aes(size=Frequency,col=Model)) + geom_text_repel(data=df[(df["adjusted_pvalue"]<p_value_threshold) ,], aes(label=names,fontface=8,size=8))+scale_size_continuous(range = c(3,12)) + 
             scale_color_manual(values=colors_analysis,limits=names(colors_analysis))+theme_bw()+  ## respect color in feature importance
            theme(plot.title = element_text(hjust = 0.5,size=35),axis.text=element_text(size=22),axis.title=element_text(size=size_title,face="bold"))  + theme(legend.position=legend.position,legend.text=element_text(size=25),legend.title=element_text(size=25)) + 
#scale_x_continuous(limits = quantile(df$beta, quantile,na.rm = T))+
    ylab(ylab)+xlim(-max(abs(df$beta),na.rm=T),max(abs(df$beta),na.rm=T))+geom_hline(yintercept=-log(p_value_threshold,10),linetype="dashed")+
    geom_vline(xintercept=0,linetype="dashed")+
    xlab(expression(beta))+ggtitle(target)
    return(s)

}
              
              
multivariate_ranking <- function (df_final, targets , covariates,num_iterations=10) {
    registerDoMC(cores=50)
    set.seed(9)
    data <- df_final
    df <- NULL
    for (var in targets){
        df_multi <- NULL
        i <- 1
        for (i in c(1:num_iterations)){
            if(var=="hb"){
                data <- data[data$hb<20,]
            }
                set.seed(i)
                res1 <- cv.glmnet(data.matrix(data[,covariates]), data[,var], 
                                  family=ifelse(var %in% c("gender","ahd"),"binomial","gaussian")
                                  ,alpha=1,nfolds=10,parallel =T)
                df_multi <- cbind(df_multi,as.matrix(coef(res1,s="lambda.min")))
                i <- i+1
                }
        l <- data.frame(coef = rowSums(df_multi)/num_iterations)
        l$names <- str_remove(rownames(l),"full_component_")
        l <- l[-1,]
        l$abs_coef <- abs(l$coef)
        l$coef <- l$coef/sum(l$abs_coef)
        l$abs_coef <- l$abs_coef/sum(l$abs_coef)
        l <- l[order(l$abs_coef,decreasing=T),]
        l$rank <- 1: nrow(l)
        l[l$coef==0,"rank"] <- "-"
        l$var <- var
        df <- rbind(df,l)
        df$names <- factor(df$names,levels=str_remove(covariates,"full_component_"))
    }
    ggplot(df, aes(factor(names),var, fill = coef))+
        geom_tile(color = "white")+
         scale_fill_gradientn(colours = c("firebrick3", "white", "dodgerblue4"),
                       limits = c(-max(l$abs_coef), max(l$abs_coef)),
        oob = squish)+
        coord_fixed()+geom_text(aes(names,var, label = rank), color = "black", size =5) +
        xlab("")+
        ylab("")+
        theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1,face="bold"),
        axis.text.y = element_text( vjust = 1, size = 12, hjust = 1,face="bold"), panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.justification = c(2, 0),
        legend.direction = "horizontal",
        legend.text=element_text(size=6),
        legend.title=element_text(size=15),
        plot.title=element_text(size=20,face="bold",hjust=0.5))+
        guides(fill = guide_colorbar(barwidth = 5, barheight = 4,title.position = "top", title.hjust = 0.5),
        legend.title=element_text(size=15,face="bold"))+
        ggtitle(paste("Multivariate coefficients Lasso ranking repeated ",sep= as.character(num_iterations), " times "))+
        labs(fill="Relative coefficients")
        }
              
              
multivariate_regression <- function(df_final,target,features=c(gen,cyto),fam="gaussian",num_iterations=100,lambda=NULL,top_n=10,size_title=25,size_axis=12,
                                    colors_analysis=colors_analysis){                  
    registerDoMC(cores=50)


    data <- df_final
    df_multi <- NULL
    i <- 1
    for (i in c(1:num_iterations)){
        res1 <- cv.glmnet(data.matrix(data[,features]), data[,target], family=fam,alpha=1,nfolds=10,parallel =T , lambda=lambda,seed=i)
        df_multi <- cbind(df_multi,as.matrix(coef(res1,s="lambda.min")))
        i <- i+1
        }
    l <- data.frame(coef = rowSums(df_multi)/num_iterations)
    l$names <- rownames(l)
    l <- l[-1,]
    l$abs_coef <- abs(l$coef)
    l$Model <- ifelse(l$names %in% gen,"gen",
                  ifelse(l$names %in% cyto,"cyto",
                        ifelse(l$names %in% clin, "clin",
                              ifelse(l$names %in% demo, "demo",
                                    ifelse(l$names %in% eln, "eln",
                                           ifelse(l$names %in% mrd,"mrd","comp"))))))
    l <- l %>% top_n(top_n,abs_coef)
    t <- ggplot(l[l$coef!=0,],aes(x=reorder(names,abs_coef),y=coef,fill=Model))+geom_bar(stat="identity")+coord_flip()+theme_bw()+
    theme(plot.title = element_text(hjust = 0.5,size=25),axis.text=element_text(size=size_axis),axis.title=element_text(size=size_title,face="bold"),
         legend.position="none")+
    scale_fill_manual(values=colors_analysis)+ylab(paste(paste("Lasso Regression Coefficients for ",target,sep="")," (averaged over 100 repetitions)",sep="")) + xlab("Top 10 Features")
    #     png(paste(target,"_lasso.png",sep=""),width=4500,height=2500,res=200)
    return(t)
    #     dev.off()
}
                    

####################################### Heatmap ####################################### 
heatmaps <- function(df_final, type=c("not_ordered","ordered"),cols_to_plot=c(genes,cytos),
                    col_order=c("t(15;17)","inv(16)","t(8;21)","t(v;11)","t(9;11)","t(6;9)","inv(3)","CEBPA_bi","CEBPA_mono","NPM1",
                                "TP53","complex","-5","-7","-17","-9","-13","-20","-18","-16","-12","-3","-Y","+8","+13","+21","+22",
                                "+11","ZRSR2", "U2AF1_p.S34","U2AF1_p.Q157", "SRSF2", "SF3B1", "SF1", "NF1", "CUX1","ASXL1", "STAG2",
                                "BCOR", "MLL", "EZH2","PHF6","RUNX1","SETBP1","DNMT3A","IDH1","IDH2_p.R140","IDH2_p.R172","WT1","ITD",
                                "TET2","FLT3_TKD","FLT3_other")){
    
    cols_to_keep <- c(cols_to_plot,"principal_component_numeric")
    features_cols <- colnames(df_final[,cols_to_plot])
    tmp_all <-df_final[,cols_to_keep]
    for (ord in type){
    #     png(paste(ord,"_heatmap.png"),width=3000,height=3500,res=200)
        if (ord=="ordered"){
            tmp_all <- df_final[,cols_to_keep][order(df_final[,cols_to_keep]$principal_component_numeric),]
        }

        transp_df_all <-(as.data.frame(t(tmp_all)))
        transp_df_all$features <- rownames(transp_df_all)
        transp_df_all <- transp_df_all[! row.names(transp_df_all) == "principal_component_numeric",]
        col_order <- c(col_order,setdiff(features_cols,col_order))
        col_order <- rev(col_order)
        transp_df_all <- transp_df_all[col_order,]
        melt.data_all<-melt(setDT(transp_df_all),id.vars ="features", variable_name=colnames(transp_df_all))
        melt.data_all$features <- factor(melt.data_all$features,levels=col_order)  # rev because otherwise the first are at the bottom


        pal <- ifelse(is.element(melt.data_all$features,genes),"#2b8cbe","#BFBFBF") 
        col_sep <- tabulate(df_final[,cols_to_keep]$principal_component_numeric)   ## to separate by component
        if(ord=="ordered"){
            cat("Order of Components :","\n")
            cat(paste(unique(df_final[order(df_final[,cols_to_keep]$principal_component_numeric),]$principal_component)," ,"))
        }
            # plot the heatmap
        p <- ggplot(melt.data_all, aes(variable,features)) +  
        geom_raster(aes(fill = factor(value)), show.legend = FALSE)+
        scale_fill_manual(values = c("0" = "grey45", "1" = "#e79f00"))+
        theme(axis.text.y = element_text(colour = pal,size=20,face="bold"))+
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
        ylab("")
        if(ord=="ordered"){
            p <- p+geom_vline(xintercept = cumsum(col_sep)[c(1,2,3,4,6,7,8,9,10)] + 0.5, col = "#009E73", linetype = 1, size = 1) +
            geom_vline(xintercept = cumsum(col_sep)[c(5,11,12,13,14)] + 0.5, col = "#CC79A7", linetype = 2, size = 1) +
            geom_vline(xintercept = cumsum(col_sep)[c(15,16)] + 0.5, col = "white", linetype = 3, size = 1) 
        } else {
            p <- p+theme(axis.text.y=element_blank())
        }
        
        return(p) 
        }
            
} 
              
####################################### Genetic landscape of strats ####################################### 
              
              
stratas_features_frequency <- function(df_final="",genes="",cytos="",stratas="",cols_to_keep=c(gen,cyto,"inv(3)","t(6;9)","t(9;11)","t(15;17)"),facet_type=T,colors="",
                                   text.y=20){
  

    set_notebook_plot_size(60,10*length(stratas))

    columns <- c(cols_to_keep,stratas)

    df_final_selected_columns <- df_final[,columns]
    stratas_repartition <- data.frame(category = colnames(df_final_selected_columns))


    # Get frequency of genes and cytos for each selected columns
    for (strata in stratas){

        stratas_repartition[sprintf('%s', strata)] <- apply(stratas_repartition, 1, 
                                                            function(s) (sum(df_final_selected_columns[df_final_selected_columns[,strata] == 1,
                                                                         s['category']])/nrow(df_final_selected_columns[df_final_selected_columns[,strata] ==1,])))
    }

    # Remove stratas rows
    stratas_repartition <- stratas_repartition[!stratas_repartition$category %in% stratas,]


    #Transpose stratas_repartition and put category row as column and remove it
    transpose_stratas_repartition <- t(stratas_repartition)
    colnames(transpose_stratas_repartition) = transpose_stratas_repartition["category", ]                                                        
    transpose_stratas_repartition <- transpose_stratas_repartition[-1,]

    names <- rownames(transpose_stratas_repartition)
    rownames(transpose_stratas_repartition) <- NULL
    data <- as.data.frame(cbind(names,transpose_stratas_repartition))



    data$names <- str_remove(data$names,"principal_component_")
    data$names <- ifelse(data$names=="TP53_complex","TP53-complex",
                                            ifelse(data$names=="DNMT3A_IDH1_2","DNMT3A_IDH",data$names))
    data <- data[!(data$names=="no_events"),]

    melt_data <- melt(setDT(data),id="names",name="Frequency", variable.name="Variable")
    melt_data$value <- as.numeric(melt_data$value)  
    melt_data$fact <- factor(melt_data$names, levels=c(data$names))
    melt_data$separation <- ifelse(melt_data$Variable %in% genes,"genes",
                                    ifelse(melt_data$Variable %in% c("t(15;17)","t(v;11)","t(9;11)",'t(8;21)','t(6;9)','t(8;21)',"others_transloc","inv(3)","inv(16)"),"fusions","CNAs"))
    melt_data$separation <- factor(melt_data$separation, levels=c("genes","CNAs","fusions"))

    p <- ggplot(melt_data, aes(x=Variable, y=value, fill=fact))+
      geom_bar(stat="identity")
    if(facet_type==T){
        p <- p + facet_grid(fact~separation,scales = "free_x", space = "free_x")
    }else{
        p <- p + facet_grid(fact~.)
    }                                                           

    p+theme_minimal()+theme(strip.text.x = element_text(face="bold", size=34,lineheight=5.0),strip.text.y = element_blank(),
                            #strip.text.y = element_text(face="bold",size=18,lineheight=5.0),
    strip.background = element_rect( colour="black",size=1),legend.position="none")+
        ylab("")+xlab("")+theme(plot.title = element_text(hjust = 0.5,size=120,face="bold"),axis.text.x = element_text(angle = 90, hjust = 1, size=ifelse(length(stratas)>1,40,40),face="bold",vjust=0.5),axis.text.y = element_text(size=text.y))+
         theme(axis.title.y = element_text(size = 40, angle = 90, vjust = 5,face="bold"))+
            scale_fill_manual(values=colors)+theme(legend.position = "none")+
            scale_y_continuous(breaks=seq(0, 1, 0.2),position="right")

}
              
eln_stratas_freq <- function(df_final="",cols="",flip=T){     
    p <-list()
    for (co in cols){
        tmp <- df_final[df_final[,co]==1,]
        tmp$eln <- ifelse(tmp$eln_2017==1,"adverse",ifelse(tmp$eln_2017==2,"intermediate","favorable"))
        p[[co]] <- (ggplot(tmp,aes(x="")) +
            geom_bar(aes(y=..count../sum(..count..),fill=factor(eln)))+scale_fill_manual(values=c(adverse="#EE9937",favorable="#2b8cbe",intermediate="#59A08B"))+
            theme_minimal()+theme(axis.text=element_text(size=12),legend.position="none",axis.title=element_text(size=14,face="bold"))+ylab("")+xlab(""))+
        scale_y_continuous(breaks=seq(0, 1, 0.5),position="right")
        if(flip==T){
            p[[co]] <- p[[co]] + coord_flip()
        }
        
    }
    return (p)
}                
                                                                                                                       
####################################### End Genetic landscape of strats #######################################                                                            

                                                            
####################################### Feature Importance #######################################                                                            
                                                            
genes_cytos <- c(genes,cytos)
clin_demo <- c(clin,demo)
pal <- c(eln="#EE9937",comp="#59A08B",genes_cytos="#BFBFBF",clin_demo="#870C14")
importances <- function(combinations,top_n="all",remove_spec_gen_cyto=F,all_algo=T,algo="",path="Prognosis/Importance/"){
    p <- list()
    for (comb in combinations){
        df_reshuffling <- read.table(paste(path,"_reshuffle_importance.tsv",sep=comb))
        df_reshuffling$Ratio <- df_reshuffling$ref_CI/df_reshuffling$permuted_CI
        df_reshuffling <- df_reshuffling[,c("feature","Ratio","algo","model")]
        df_reshuffling$feature <- as.character(df_reshuffling$feature)
        df_reshuffling[grep("^add_",df_reshuffling$feature),"feature"] <- str_replace(df_reshuffling[grep("^add_",df_reshuffling$feature),"feature"],"add_","+")
        df_reshuffling[grep("^del_",df_reshuffling$feature),"feature"] <- str_replace(df_reshuffling[grep("^del_",df_reshuffling$feature),"feature"],"del_","-")
        df_reshuffling[df_reshuffling$feature=="minusy","feature"] <- "-Y"
        df_reshuffling[grep("^t_",df_reshuffling$feature),"feature"] <- paste(str_replace(str_replace
                                                                        (df_reshuffling[grep("^t_",df_reshuffling$feature),"feature"],"t_","t("),"_",";"),")",sep="")
        df_reshuffling[grep("^inv_",df_reshuffling$feature),"feature"] <- paste(str_replace
                                                                        (df_reshuffling[grep("^inv_",df_reshuffling$feature),"feature"],"inv_","inv("),")",sep="")
        

        df_reshuffling$feature <- ifelse(df_reshuffling$feature=="principal_component_t_8_21","principal_component_t(8;21)",
                                        ifelse(df_reshuffling$feature=="principal_component_t_6_9","principal_component_t(6;9)",
                                            ifelse(df_reshuffling$feature=="principal_component_t_15_17","principal_component_t(15;17)",
                                                ifelse(df_reshuffling$feature=="principal_component_t_11","principal_component_t(11)",
                                                       ifelse(df_reshuffling$feature=="principal_component_inv_3","principal_component_inv(3)",
                                                              ifelse(df_reshuffling$feature=="principal_component_inv_16","principal_component_inv(16)",
                                                                     ifelse(df_reshuffling$feature=="principal_component_CEBPA_bi","principal_component_CEBPAbi",
                                                                  df_reshuffling$feature)))))))

        
        if(remove_spec_gen_cyto){  #useful for redundancy class defining genes and cytos vs those classes
            df_reshuffling <-df_reshuffling[!df_reshuffling$feature %in% c("NPM1","t(8;21)","inv(3)","inv(16)","t(6;9)","WT1","t(15;17)"),]
        }
        
        df_reshuffling <- df_reshuffling[,c('feature','Ratio','algo')] %>% group_by(feature,algo)
        df_reshuffling$feature <- str_replace(df_reshuffling$feature,"principal_component_","class ")
        df_reshuffling$algo <- ifelse(df_reshuffling$algo=="RFX","Random Effects",
                                      ifelse(df_reshuffling$algo=="RFS","Random Forest",
                                             ifelse(df_reshuffling$algo=="Ridge","Ridge",
                                                    ifelse(df_reshuffling$algo=="Lasso","Lasso",
                                                           ifelse(df_reshuffling$algo=="Elastic_net","Elastic-net","Cox Vanilla")))))
        
        df_reshuffling$algo <- factor(df_reshuffling$algo, levels = c("Random Forest","Random Effects","Lasso","Ridge","Elastic-net","Cox Vanilla"))
        
        df_reshuffling$category <- ifelse(is.element(df_reshuffling$feature,eln),"eln",
                                               ifelse((grepl("class",df_reshuffling$feature) | df_reshuffling$feature=="overlap"),"comp",
                                                                 ifelse(is.element(df_reshuffling$feature,clin_demo),"clin_demo","genes_cytos")))
        name_comb <-ifelse(comb=="eln","ELN 2017",
                         ifelse(comb=="comp","Classes",
                                ifelse(comb=="comp_ITD","Classes-ITD",
                                       ifelse(comb=="comp_all_gen_cyto","Classes+Genes+Cytos",
                                       ifelse(comb=="comp_ITD_gen_cyto","Classes-ITD+Genes+Cytos",       
                                       ifelse(comb=="gen_cyto","Genes+Cytos",
                                       ifelse(comb=="eln_all_gen_cyto_clin_demo","ELN+Genes+Cytos+Clin+Demo",
                                              ifelse(comb=="comp_clin_demo","Classes+Clin+Demo",
                                                     ifelse(comb=="comp_ITD_clin_demo","Classes-ITD+Clin+Demo",
                                                    ifelse(comb=="clin_demo","Clin+Demo",
                                                           ifelse(comb=="all_gen_cyto","Genes+Cytos",
                                                                  ifelse(comb=="gen" | comb=="all_gen","Genes",
                                                                         ifelse(comb=="eln_all_gen_cyto","ELN+Genes+Cytos",
                                                                                ifelse(comb=="eln_all_gen","ELN+Genes",
                                                                         ifelse(comb=="cyto" | comb=="all_cyto","Cytos",comb)))))))))))))))
        if(all_algo==T){
            df_reshuffling <- df_reshuffling[df_reshuffling$algo %in% c("Random Forest","Random Effects","Lasso","Ridge","Elastic-net"),]
        } else {
            df_reshuffling <- df_reshuffling[df_reshuffling$algo %in% algo,]
        }
        
        title <- paste("Feature importance for",name_comb,sep=" ")
        a <- group_by(df_reshuffling, feature) %>% summarize(m = max(Ratio))
        if(top_n!="all"){
            df_reshuffling <- df_reshuffling[df_reshuffling$feature %in% a[order(a$m,decreasing=T),]$feature[1:top_n],]
            title <- paste("Top",top_n,sep=" ",paste(" Feature importance for",name_comb,sep=" "))
        }
        df_reshuffling$feature <- factor(df_reshuffling$feature,levels=a[order(a$m),]$feature)
        #
        if(comb=="all_gen_cyto_comp"){
            df_reshuffling <- df_reshuffling[!df_reshuffling$feature %in% c("t(15;17)","inv(16)","t(8;21)","t(6;9)","inv(3)","NPM1","CEBPA_bi","WT1"),]
        }

        p[[comb]] <- ggplot(na.omit(df_reshuffling),aes(x=feature,y=Ratio,fill=category))+theme_bw()+
        geom_boxplot()+facet_wrap(~algo,ncol=6,scales="free_x")+theme(plot.title = element_text( size=50, face="bold",hjust=0.5,vjust=1.5),
                                                                      legend.text=element_text(size=20),legend.position="top",
                  axis.text.x = element_text(angle = 90, hjust = 0.8,vjust=0.8,size=28),
                                                                      axis.text.y=element_text(size=28),
                  axis.title = element_text(size = 35),strip.text = element_text(size=35))+xlab("Selected Features") +
                ylab("Ratio (ref_CI/Permuted_CI) ") + ggtitle(title) + theme(legend.position = "none") +
                scale_fill_manual(values=pal,limits=names(pal))
        p[[comb]] <- p[[comb]] + coord_flip()
        
    }
    return(p)
}
                                                            
####################################### End Feature Importance #######################################
                                                            
####################################### Model Comparison #######################################
                                                            
comparison <- function(path_training="Prognosis/Results/", models = c("eln","comp"), remove_RFX_ = F ,tests="",color="",FUN=median,compare_algo=F,label.y=0.75,add_val=FALSE,val_only=FALSE,add_colors=T ){

    with_colors = T

    df_comparison <- c()
    for (mod in models){
        
        tmp <- cbind(read.table(paste(paste(path_training,mod,sep=""),".tsv",sep="")),SVM=SVM[,mod])
        tmp$model <- ifelse(mod=="eln","ELN 2017",
                 ifelse(mod=="comp","Classes",
                 ifelse(mod=="comp_ITD","Classes-ITD",
                 ifelse(mod=="eln_gen_cyto_clin_demo","ELN+Genes+Cytos+Clin+Demo",
                 ifelse(mod=="comp_clin_demo","Classes+Clin+Demo",
                 ifelse(mod=="comp_ITD_clin_demo","Classes-ITD+Clin+Demo",
                 ifelse(mod=="clin_demo","Clin+Demo",
                 ifelse(mod=="all_gen_cyto","Genes+Cytos",
                 ifelse(mod=="gen_cyto","Genes+Cytos",
                 ifelse(mod=="gen","Genes",
                 ifelse(mod=="all_gen","Genes",
                 ifelse(mod=="cyto","Cytos",
                 ifelse(mod=="eln_cyto","ELN+Cytos",
                 ifelse(mod=="eln_gen","ELN+Genes",
                 ifelse(mod=="eln_gen_cyto","ELN+Genes+Cytos",
                 ifelse(mod=="gen_cyto_clin_demo","Genes+Cytos+Clin+Demo",
                 ifelse(mod=="all_gen_cyto_clin_demo","Genes+Cytos+Clin+Demo",mod)))))))))))))))))
        
        if(with_colors){
                with_colors = is.element(unique(tmp$model),names(colors))
        }
        tmp$Cohort <- "MRC"
        
        if (add_val==T) {
                tmp_val <-read.table(paste(paste("Validation/Prognosis/",mod,sep=""),".tsv",sep=""))
                tmp_val$Cohort <- "Validation"
            tmp_val$model <- ifelse(mod=="eln","ELN 2017",
                 ifelse(mod=="comp","Classes",
                 ifelse(mod=="comp_ITD","Classes-ITD",
                 ifelse(mod=="eln_gen_cyto_clin_demo","ELN+Genes+Cytos+Clin+Demo",
                 ifelse(mod=="comp_clin_demo","Components+Clin+Demo",
                 ifelse(mod=="comp_ITD_clin_demo","Classes-ITD+Clin+Demo",
                 ifelse(mod=="clin_demo","Clin+Demo",
                 ifelse(mod=="all_gen_cyto","Genes+Cytos",
                 ifelse(mod=="gen","Genes",
                 ifelse(mod=="all_gen","Genes",
                 ifelse(mod=="cyto","Cytos",
                 ifelse(mod=="eln_cyto","ELN+Cytos",
                 ifelse(mod=="eln_gen","ELN+Genes",
                 ifelse(mod=="eln_gen_cyto","ELN+Genes+Cytos",
                 ifelse(mod=="gen_cyto_clin_demo","Genes+Cytos+Clin+Demo",
                 ifelse(mod=="all_gen_cyto_clin_demo","Genes+Cytos+Clin+Demo",mod))))))))))))))))
                df_comparison <- rbind(df_comparison,tmp,tmp_val)
        } else if(val_only==F) {
                df_comparison <- rbind(df_comparison,tmp)
        } else {
            tmp_val <-read.table(paste(paste("Validation/Prognosis/",mod,sep=""),".tsv",sep=""))
            tmp_val$Cohort <- "Validation"
            tmp_val$model <- ifelse(mod=="eln","ELN 2017",
                 ifelse(mod=="comp","Classes",
                 ifelse(mod=="comp_ITD","Classes-ITD",
                 ifelse(mod=="eln_gen_cyto_clin_demo","ELN+Genes+Cytos+Clin+Demo",
                 ifelse(mod=="comp_clin_demo","Classes+Clin+Demo",
                 ifelse(mod=="comp_ITD_clin_demo","Components-ITD+Clin+Demo",
                 ifelse(mod=="clin_demo","Clin+Demo",
                 ifelse(mod=="all_gen_cyto","Genes+Cytos",
                 ifelse(mod=="gen","Genes",
                 ifelse(mod=="all_gen","Genes",
                 ifelse(mod=="cyto","Cytos",
                 ifelse(mod=="eln_cyto","ELN+Cytos",
                 ifelse(mod=="eln_gen","ELN+Genes",
                 ifelse(mod=="eln_gen_cyto","ELN+Genes+Cytos",
                 ifelse(mod=="gen_cyto_clin_demo","Genes+Cytos+Clin+Demo",
                 ifelse(mod=="all_gen_cyto_clin_demo","Genes+Cytos+Clin+Demo",mod))))))))))))))))
            df_comparison <- rbind(df_comparison,tmp_val)
        }
        
        
        
        
    }
    df_comparison <- melt(df_comparison,id.vars=c('model','Cohort'))
    if (remove_RFX_==T){df_comparison <- df_comparison[df_comparison$variable!= "RFX_",]}
    df_comparison $title <- "Detailed Comparison"
    df_comparison$title_bis <- "Feature Combinations"
    df_comparison$title3 <- "Algorithms Comparison on top performing feature combination"
    
    lis <- c("CoxGLM_0","CoxGLM_1","CoxGLM_0.4","CoxGLM_0.8","CoxBoost","RFS_1050_20","SVM")
    p1 <- ggplot(df_comparison[df_comparison$variable %in% lis,],aes(x=variable,y=value)) + geom_boxplot(aes(fill = model), notch=T)+theme_bw()+
    theme(strip.text = element_text(size=30),legend.position="top",axis.text.x = element_text(angle = 35,
    hjust = 0.8,vjust=0.8,size=20),axis.text.y = element_text(size = 20),
         axis.title.x=element_text(size=30),axis.title.y=element_text(size=30),legend.title=element_text(size=30),legend.text=element_text(size=20)) +
    ylab("C-Index") + xlab("") + facet_wrap( ~ title) + 
        scale_fill_manual(values=c("olivedrab3","darkorange","thistle4","lightsalmon","lightsteelblue","cornflowerblue","lightcyan3","paleturquoise3"))
    
    p3 <- ggplot(df_comparison[df_comparison$variable %in% lis,],aes(x=reorder(variable,value),y=value)) + geom_boxplot(aes(fill = model), notch=T)+theme_bw()+
    theme(strip.text = element_text(size=30),legend.position="top",axis.text.x = element_text(angle = 25, hjust = 0.8,vjust=0.8,size=20),
          axis.text.y = element_text(size = 18),
         axis.title.x=element_text(size=30),axis.title.y=element_text(size=25),
          legend.title=element_text(size=30),legend.text=element_text(size=18)) +
    ylab("C-Index") + xlab("") + facet_wrap( ~ title3) + 
        scale_fill_manual(values=c("lightsteelblue"))
                          
                          
    p2 <- ggplot(df_comparison[df_comparison$variable %in% lis,], aes(x=reorder(model,value,FUN=FUN),y=value,fill=Cohort))+theme_bw()
    if(add_colors==T){
        p2 <- p2 + geom_rect(fill = 'lightsteelblue', xmin = 0, xmax = 5.5, ymin =0, ymax = Inf, alpha =1)+
        geom_rect(fill = 'lavenderblush', xmin = 5.5, xmax = Inf, ymin =0, ymax = Inf, alpha =1)
    } 
    p2 <- p2 + geom_boxplot(fatten=3)+scale_fill_manual(values=c("#708090"))+
    theme(strip.text = element_text(size=25),axis.text.x = element_text(angle = 35,
    hjust = 0.8,vjust=0.8,size=15),axis.text.y = element_text(size = 15),
             axis.title.x=element_text(size=30),axis.title.y=element_text(size=30),legend.title=element_text(size=30),
          legend.text=element_text(size=20))+
     stat_compare_means(paired=F,comparison=tests,method="wilcox.test",label="p",size=5,label.y=label.y,color="blue")+
    xlab("")+ylab("C-Index")+ facet_wrap( ~ title_bis) 
    if(compare_algo==F){
        return(list(p1,p2))
    }else{
        return(p3)
    }
}
                                                            
                                                            
                                                            
                                                            
                                                            
####################################### Transplant Preparation #######################################
                                                            
prepare_transplant_data <- function(data,columns){
    data$ID <- 1:nrow(data)
    columns_matrix <- as.matrix(data[,columns])
    tdc_No_TPL <- NULL 
    tdc_TPL <- NULL 

    count <- 0
    for(pat_ID in 1:nrow(data)){
        
     # Transplanted patients
        if(!is.na(data$TPL_date[pat_ID])){  
        tdc_No_TPL <- rbind(tdc_No_TPL, 
                      c(data$ID[pat_ID], 0,data$TPL_date[pat_ID]-data$CR[pat_ID] , 0,0, columns_matrix[pat_ID,]))

        tdc_TPL <- rbind(tdc_TPL, 
                      c(data$ID[pat_ID], data$TPL_date[pat_ID]-data$CR[pat_ID],data$OS[pat_ID]-data$CR[pat_ID] , data$OS_stat[pat_ID],1, columns_matrix[pat_ID,]))

     # No Transplant    
        }else{
        tdc_No_TPL <- rbind(tdc_No_TPL, 
                      c(data$ID[pat_ID], 0,data$OS[pat_ID]-data$CR[pat_ID] , data$OS_stat[pat_ID],0,  columns_matrix[pat_ID,]))
        }
    }

    tdc_No_TPL <- data.frame(tdc_No_TPL)
    tdc_TPL <- data.frame(tdc_TPL)

    colnames(tdc_No_TPL) <- colnames(tdc_TPL) <- c("ID", "Start", "Stop", "Status", "Transplant", columns)
    
    return (list(No_TPL=tdc_No_TPL,TPL=tdc_TPL))
}


                  
prepare_transplant_pvalues <- function(data,columns,column_tdc){
    data$ID <- 1:nrow(data)
    columns_matrix <- as.matrix(data[,columns])
    
    tdc_No_TPL <- NULL 
    tdc_TPL <- NULL 

    count <- 0
    for(pat_ID in 1:nrow(data)){
        
     # Transplanted patients
        if(!is.na(data$TPL_date[pat_ID])){  
        tdc_No_TPL <- rbind(tdc_No_TPL, 
                      c(data$ID[pat_ID], 0,data$TPL_date[pat_ID]-data$CR[pat_ID] , 0,0, columns_matrix[pat_ID,],data[,column_tdc][pat_ID]+data$CR[pat_ID]))

        tdc_TPL <- rbind(tdc_TPL, 
                   c(data$ID[pat_ID], data$TPL_date[pat_ID]-data$CR[pat_ID],data$OS[pat_ID]-data$CR[pat_ID] ,data$OS_stat[pat_ID],1, columns_matrix[pat_ID,],
          data[,column_tdc][pat_ID]+data$TPL_date[pat_ID]))

     # No Transplant    
        }else{
        tdc_No_TPL <- rbind(tdc_No_TPL, 
                      c(data$ID[pat_ID], 0,data$OS[pat_ID]-data$CR[pat_ID] , data$OS_stat[pat_ID], 0,  columns_matrix[pat_ID,],data[,column_tdc][pat_ID]+data$CR[pat_ID]))
        }
    }

    tdc_No_TPL <- data.frame(tdc_No_TPL)
    tdc_TPL <- data.frame(tdc_TPL)

    colnames(tdc_No_TPL) <- colnames(tdc_TPL) <- c("ID", "Start", "Stop", "OS_Status", "Transplant", columns,paste0(column_tdc,"_tdc"))
    
    return (list(No_TPL=tdc_No_TPL,TPL=tdc_TPL))
}

####################################### End Transplant Preparation #######################################


####################################### Helper for calculator ########################################### 

prep_probtrans_df <- function(obj,
                              from,
                              ord,
                              state_names,
                              conf.type = "log",
                              conf.int = 0.95) {
  
  # For R cmd check
  . <- state_num <- state <- se_state <- prob <- se <- cum_probs <-  NULL
  
  # Read in probtrans object
  df <- data.table::data.table(obj[[from]])
  
  # Subset states with at least one non zero probability
  zero_prob_cols <- apply(df, 2, function(col) !all(col == 0))
  
  # Condition pstate cols
  pstate_cols <- grepl(x = names(df), pattern = "pstate|time")
  condition_pstate <- zero_prob_cols & pstate_cols  

  # Prepare 
  df_pstate <- data.table::melt.data.table(
    data = df[time != Inf, .SD, .SDcols = condition_pstate],
    id.vars = "time",
    variable.name = "state",
    value.name = "prob"
  ) 
  df_pstate[, state_num := as.numeric(gsub(x = state, pattern = "pstate", replacement = ""))]
  df_long <- df_pstate
  #df_long[, ':=' (CI_low = NA, CI_upp = NA)]
  
  # Check if standard errors were computed (only time)
  se_cols <- grepl(x = names(df), pattern = "se|time")
  if (sum(se_cols) > 1) {
    condition_se <- zero_prob_cols & se_cols
    
    # Se's long
    df_se <- data.table::melt.data.table(
      data = df[time != Inf, .SD, .SDcols = condition_se],
      id.vars = "time",
      variable.name = "se_state",
      value.name = "se"
    ) 
    df_se[, state_num := as.numeric(gsub(x = se_state, pattern = "se", replacement = ""))]
    
    # Put se's and pstates together
    df_long <- data.table::merge.data.table(
      x = df_pstate, 
      y = df_se, 
      by = c("time", "state_num")
    )
    
    # Add CI for probabilities, do this on log
    df_long[, ':=' (
      CI_low = make_prob_confint(prob, se, conf.type, conf.int, bound = "low"),
      CI_upp = make_prob_confint(prob, se, conf.type, conf.int, bound = "upp")
    )]  #, by = .(time, state)] %>% 
  }
  
  # Order and label factor with state names
  df_long[, state := factor(
    x = state_num, 
    levels = ord, 
    labels = state_names[ord]
  )]
  
  data.table::setorder(df_long, time, state)
  df_long[, cum_probs := cumsum(prob), by = time]
  
  # Compute upper and lower bounds of ribbons
  df_long[, ':=' (
    low = c(0, cum_probs[-length(cum_probs)]),
    upp = c(cum_probs[-length(cum_probs)], 1)
  ), by = time] 
  
  # Df with steps for ribbon
  df_steps <- rbind(df_long, df_long) 
  data.table::setorder(df_steps, time)
    
  # Shift time by 1, creating steps, equi to dplyr::lead(time, n = 1)
  df_steps[, time := data.table::shift(
    x = time, 
    fill = NA, 
    n = 1, 
    type = "lead"
  ), by = state] 
 
  return(df_steps[!is.na(time)]) 
}


make_labelled_plot <- function(df_steps,
                               state_names,
                               xlim,
                               cex,
                               lwd) {
  
  # For R cmd check
  . <- state <- low <- upp <- mean_cprob <- label <- NULL
  
  # Max follow-up, depends on xlim
  time_vec <- df_steps[df_steps$time <= xlim[2], ]$time
  max_t <- time_vec[length(time_vec)]

  # Prep labels - at end of follow-up
  df_labels <- df_steps[, ':=' (
    label = ifelse(time == max_t, as.character(state), NA_character_),
    mean_cprob = (low + upp) / 2
  )]  
  
  # Begin plots
  p <- ggplot2::ggplot(
    data = df_labels[time <= max_t, .SD[-.N], by = state],
    ggplot2::aes(time, mean_cprob)
  ) + 
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = low, ymax = upp, fill = state),
      col = "black", size = lwd
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = label, x = time), 
      na.rm = T,
      hjust = 1, 
      size = cex
    ) 
    
  return(p)
}


make_prob_confint <- function(prob, 
                              se, 
                              conf.type = "log",
                              conf.int = 0.95, 
                              bound) {
  # Get critical value
  if (!is.null(conf.int)) {
    crit <- qnorm((1 - conf.int) / 2, lower.tail = FALSE)
  }

  if (conf.type == "log") {
    low <- exp(log(prob) - crit * se / prob)
    upp <-  exp(log(prob) + crit * se / prob)
      
  } else if (conf.type == "plain") {
    low <- prob - crit * se
    upp <- prob + crit * se
    
  } else {
    low <- prob
    upp <- prob
  }
  
  # Check no vals above 1 or below zero
  low <- ifelse(low < 0, 0, low)
  upp <- ifelse(upp > 1, 1, upp)
  
  # Return upper of lower bound
  if (bound == "upp") {
    return(upp) 
  } else return(low)
}
                          
                          
plot_probtrans_Conf_Interval <- function(probtrans_patient,colors,ord,
                                         from=1,xlim=c(0,10.5),ylim=c(0, 1),conf.int=0.95,conf.type="log",col_ribb="grey70",
                                        xlab="Time",ylab="Probability",main="",state_plot=F){
    x <- probtrans_patient

    low <- upp <- state <- prob <- CI_low <- CI_upp <- NULL

    # Get names of states and state numbers

    state_names <- dimnames(x$trans)[[1]]


    # Check stacking order
    ord <- 1:length(state_names)
    cex <- 8
    lwd <- 0.5


    # Create copy of pb object so original is NOT affected by references updating
    pb_copy <- data.table::copy(x)

    # Prepare dataframe 
    df_steps <- prep_probtrans_df(
    obj = pb_copy,
    from = from,
    ord = ord, 
    state_names = state_names,
    conf.int = conf.int,
    conf.type = conf.type
    )


    . <- state <- low <- upp <- mean_cprob <- label <- NULL

    # Max follow-up, depends on xlim
    time_vec <- df_steps[df_steps$time <= xlim[2], ]$time
    max_t <- time_vec[length(time_vec)]

    # Prep labels - at end of follow-up
    df_labels <- df_steps[, ':=' (
    label = ifelse(time == max_t, as.character(state), NA_character_),
    mean_cprob = (low + upp) / 2
    )] 

    tmp <- df_labels[time <= max_t, .SD[-.N], by = state]
    crit <- qnorm((1 - 0.95) / 2, lower.tail = FALSE)

    tmp$CI_low_bis <- exp(log(tmp$upp) - crit * tmp$se / tmp$upp)
    tmp$CI_upp_bis <- exp(log(tmp$upp) + crit * tmp$se / tmp$upp)

    tmp[tmp$CI_upp_bis>1,"CI_upp_bis"] <- 1
    tmp[tmp$CI_low_bis<0,"CI_low_bis"] <- 0

#     # Remove Conf Interval for last line (when from = 1)
# #     tmp[tmp$state=="Death without CR","CI_low_bis"] <- 10 # more than 1 so that code do not account for it
#     tmp[tmp$state=="Death without CR","CI_upp_bis"] <- 1 # more than 1 so that code do not account for it
    
#     # Remove Conf Interval for last line (when from = 2)
#     if(from==2){
#         tmp[tmp$state=="Death in Relapse","CI_low_bis"] <- 10 # more than 1 so that code do not account for it
#         tmp[tmp$state=="Death in Relapse","CI_upp_bis"] <- 10 # more than 1 so that code do not account for it
#     }
    p <- ggplot(
        data = tmp,aes(time, mean_cprob)
      ) + geom_ribbon(aes(ymin = low, ymax = upp, fill = state),
          col = "black", size = lwd)+
    geom_ribbon(aes(x=time,ymin = CI_low_bis, ymax = CI_upp_bis,group=state), 
            alpha = 0.1, 
            fill = col_ribb,
            col = col_ribb,
            na.rm = TRUE,
            linetype=2
          )+scale_fill_manual(values = colors) +
      ylab(ylab) +
      xlab(xlab) +
      coord_cartesian(expand = 0, ylim = ylim) +
      theme(legend.position = "none")+
    theme(axis.text=element_text(size=20),
            axis.title=element_text(size=20,face="bold"),
         plot.title=element_text(size=25,face="bold",hjust = 0.5))+
    scale_x_continuous(breaks=c(0,3,6,9))+
    ggtitle(main)
    
    if (state_plot){
        
        p <- p + geom_text(aes(label = label, x = time), 
          na.rm = T,
          hjust = 1, 
          size = cex)
    }
    return(p)
}
                          
                          
####################################### End Helper for calculator ########################################### 