format.pval <- function(p){
  ifelse(p<0.001,"<0.001",ifelse(p<0.05,formatC(p,3,format="f"),formatC(p,2,format="f")))
}

mySummary.onevar <- function(varname,variable,group,cont=NA,contSummary="med.90",pval.comparison=F){
  require(gdata)
  if (is.na(cont)) cont <- ifelse(is.factor(variable)|length(unique(na.omit(variable)))<=5,F,T)
  ngroup <- length(levels(group))
  mycont.summary <- function(variable,group) {
    if (is.na(match(contSummary,c("med.90","med.IQR","med.range","MIC"))))
      stop("contSummary=",contSummary," not yet implemented")
    n <- c(by(variable,group,function(x) length(na.omit(x))))
    if (contSummary=="med.90") summarystat <- by(variable,group,function(x) quantile(x,c(0.05,0.5,0.95),na.rm=T))
    if (contSummary=="med.IQR") summarystat <- by(variable,group,function(x) quantile(x,c(0.25,0.5,0.75),na.rm=T))
    if (contSummary=="med.range") summarystat <- by(variable,group,function(x) quantile(x,c(0,0.5,1),na.rm=T))
    if (!is.na(match(contSummary,c("med.90","med.IQR","med.range")))){
      summarystat.nice <- lapply(summarystat,function(x){ x <- formatC(round(x,0),0,format="f");
      paste(x[2]," (",x[1],"-",x[3],")",sep="")})
      result <- matrix("",ncol=ngroup*2+1,nrow=1)
      result[1,seq(2,ncol(result),by=2)] <- n
      result[1,seq(3,ncol(result),by=2)] <- unlist(summarystat.nice)
    }
    if (contSummary=="MIC") {
      result <- matrix("",ncol=ngroup*2+1,nrow=4)
      result[,1] <- c("","- MIC 50","- MIC 90","- range")
      result[1,seq(2,ncol(result),by=2)] <- n
      medians <- by(variable,group,function(x) quantile(x,c(0.5),na.rm=T))
      q.90s <- by(variable,group,function(x) quantile(x,c(0.9),na.rm=T))
      ranges <- by(variable,group,function(x) quantile(x,c(0,1),na.rm=T))
      result[2,seq(3,ncol(result),by=2)] <- unlist(lapply(medians,function(x){x <- formatC(round(x,0),0,format="f")}))      
      result[3,seq(3,ncol(result),by=2)] <- unlist(lapply(q.90s,function(x){x <-   formatC(round(x,0),0,format="f")}))
      result[4,seq(3,ncol(result),by=2)] <- unlist(lapply(ranges,function(x) {x <- formatC(round(x,0),0,format="f"); 
      paste(" (",x[1],"-",x[2],")",sep="")}))
      # # Replace "257" by ">256" (and "33" by ">32" for cotrimoxazole) 
      # result <- gsub("257",">256",result)
      # if (varname=="cotrimoxazole") result <- gsub("33",">32",result)
    }     
    if (pval.comparison) {
      pval <- kruskal.test(variable[group!="All patients"]~group[group!="All patients"])$p.value
      pval <- format.pval(pval)
      result <- cbind(result,"")
      result[1,ncol(result)] <- pval
    }
    result
  }
  mycat.summary <- function(variable,group) {
    ta <- table(group,variable)
    ta.prop <- ta/apply(ta,1,sum)
    ta.nice <- matrix(paste(ta,"/",apply(ta,1,sum)," (",round(100*ta.prop,0),"%",")",sep=""),nrow=nrow(ta),ncol=ncol(ta))
    result <- matrix("",ncol=ngroup*2+1,nrow=ncol(ta)+1)
    result[2:nrow(result),1] <- paste("- ",colnames(ta))
    result[2:nrow(result),seq(3,ncol(result),by=2)] <- t(ta.nice)
    result[1,seq(2,ncol(result),by=2)] <- apply(ta,1,sum) # n's
    if (pval.comparison) {
      if (ncol(table(group[group!="All patients"],variable[group!="All patients"]))==1) pval <- "NA"
      else {
        # Use simulated p-value if normal fisher-test doesn't work
        options(show.error.messages=F) 
        ft <- try(fisher.test(table(group[group!="All patients"],variable[group!="All patients"])))
        options(show.error.messages=T) 
        if (class(ft)!="try-error"){
          pval <- ft$p.value
        } else {
          warning("Simulated p-values for Fisher test used with B=10000")
          pval <- fisher.test(table(group[group!="All patients"],variable[group!="All patients"]),
                              simulate.p.value=T,B=10000)$p.value
        }  
        pval <- format.pval(pval)
      }
      result <- cbind(result,c(pval,rep("",nrow(result)-1)))
    }
    result
  }
  if (cont) r <- mycont.summary(variable,group)
  else r <- mycat.summary(variable,group)
  r[1,1] <- varname
  r
}

mySummary.allvar <- function(blvars,group,pooledGroup=F,contSummary="med.IQR",caption=NULL,filename=NA,pval.comparison=F,cont=NA){
  # contSummary can be median (90% range) "med.90" or median (IQR) "med.IQR" or median (range) "med.range" or "MIC"
  require(xtable); require(Hmisc); library(gdata)
  if (is.factor(group)) {
    levelOrder <- levels(group)
    was.factor <- T
  } else was.factor <- F
  group <- as.character(group)
  if (was.factor) group <- factor(group,levels=levelOrder[!is.na(match(levelOrder,group))])
  else group <- factor(group)
  gr.lev <- levels(group)
  if (pooledGroup){ # Add pooled summaries for all patients
    mylabels <- unlist(lapply(blvars,Hmisc::label)) # save them as rbind destroys some of them
    blvars <- rbind(blvars,blvars)
    for (i in 1:ncol(blvars)) label(blvars[,i]) <- mylabels[i] # add labels again
    group <- c(as.character(group),rep("All patients",length(group)))
    group <- factor(group,levels=c("All patients",gr.lev))
    gr.lev <- levels(group)
  }
  header1 <- c("",c(rbind(rep("",length(gr.lev)),paste(gr.lev," (N=",table(group),")",sep=""))))
  header2 <- c("Characteristic",rep(c("n","Summary statistic"),length(gr.lev)))
  result <-  rbind(header1,header2)
  if (pval.comparison) result <- cbind(result,c("Comparison"," (p-value)"))
  if (is.na(cont)[1]) cont <- rep(cont,ncol(blvars))
  for (i in 1:ncol(blvars)){
    result.i <- mySummary.onevar(varname=ifelse(Hmisc::label(blvars[,i])!="",Hmisc::label(blvars[,i]),names(blvars)[i]),
                                 blvars[,i],group,contSummary=contSummary,pval.comparison=pval.comparison,cont=cont[i])
    result <- rbind(result,result.i)
  }
  rownames(result) <- rep("",nrow(result))
  if (!is.na(filename)){  # generate html table
    x <- print(xtable(result,caption=caption),type="html",file=filename,
               caption.placement="top",include.rownames=F,include.colnames=F,
               html.table.attributes="border=1")
    cont.summary.options <- c("med.IQR","median.90","median.range","MIC")      
    cont.summary.text <- c("median (IQR)","median (90% range)","median (range)","median, 90% quantile and range")      
    x <- paste(x,"Summary statistic is absolute count (%) for categorical variables and ",
               cont.summary.text[match(contSummary,cont.summary.options)],"for continuous data.",
               ifelse(pval.comparison,"Comparisons based on Fisher's exact test for categorical data and Wilcoxon test for continuous data.",""))
    write(x,file=filename)
  }
  result
}

## Treatment comparisons overall and in subgroups - Cox regression
surv.comparison <- function(model,data,add.risk=T,add.prop.haz.test=T){
  ##---------------------------------------------------------------------------------------------------------
  ## Purpose: Summarize results for a Cox survival model with the treatment arm (variable "arm") as the main covariate
  ## !! model formula can include other covariates than arm BUT arm must be the first covariate in the model !!
  ## add.risk: if T, the event probability ("absolute risk") at time "infinity" is also displayed
  ## add.prop.haz.test: if T, a test for proportional hazards is also added
  ## Author: Marcel Wolbers, 9March2015
  ##---------------------------------------------------------------------------------------------------------
  arm.names <- levels(data$arm)
  # Table header
  header1 <- c(paste(arm.names," (n=",table(data$arm),")",sep=""),"Comparison")
  header2 <- c(rep(ifelse(add.risk,"events/n (risk [%])","events/n"),2),"HR (95%CI); p-value")
  header <- rbind(header1,header2)
  result <- rbind(header,"")
  # add number of events and risks
  fit.surv <- summary(survfit(update(model,.~arm),data=data),time=Inf,extend=T)
  events.n <- paste(fit.surv$table[,"events"],fit.surv$table[,"n.max"],sep="/")
  if (add.risk) events.n <- paste(events.n," (",formatC(100*(1-fit.surv$surv),2,format="f"),")",sep="")
  result[3,1:2] <- events.n
  # add HR, CI, p-value
  fit.coxph <- summary(coxph(model,data))
  hr <- formatC(fit.coxph$coef[1,"exp(coef)"],2,format="f")
  pval <- format.pval(fit.coxph$coef[1,"Pr(>|z|)"])
  ci <- paste(formatC(fit.coxph$conf.int[1,c("lower .95")],2,format="f"),
              formatC(fit.coxph$conf.int[1,c("upper .95")],2,format="f"),sep="-")
  hr.ci.p <- paste(hr," (",ci,"); p=",pval,sep="")
  result[3,3] <- hr.ci.p
  # add test for proportional hazards
  if (add.prop.haz.test){
    attr(model,".Environment") <- environment() # needed for cox.zph to work 
    p.prop.haz <- cox.zph(coxph(model,data))$table[1,"p"]
    result <- cbind(result,c("Test for proportional hazards","p-value",format.pval(p.prop.haz)))  
  }
  rownames(result) <- NULL
  result
}

surv.comparison.subgroup <- function(base.model,subgroup.model,data,labels,...){
  ##---------------------------------------------------------------------------------------------------------
  ## Purpose: Summarize results for a Cox survival model by treatment arm (variable "arm") and subgroup
  ## base.model: Model from which sub-group specific estimates are extracted (!! arm must be the first covariate in the model)
  ## subgroup.model: model of the form "~subgrouping.variable1+subgrouping.variable2" 
  ##                 (!! subgrouping.variable must be factors and there should be nothing on the left-hand side of the formula)
  ## ...: arguments that are passed to surv.comparison
  ## Author: Marcel Wolbers, 20March2015
  ##---------------------------------------------------------------------------------------------------------
  # result in entire population
  result <- surv.comparison(base.model,data,...)
  result <- cbind(c("Subgroup","","All patients"),result,c("Test for heterogeneity","p-value",""))
  # Preparation of models and data
  subgroup.char <- attr(terms(subgroup.model),"term.labels")
  for (k in 1:length(subgroup.char)){
    main.model <- update(base.model,as.formula(paste(".~.+",subgroup.char[k],sep=""))) 
    ia.model <- update(base.model,as.formula(paste(".~.+arm*",subgroup.char[k],sep=""))) 
    data$.subgroup.var <- data[,subgroup.char[k]]
    factor.levels <- levels(data[,subgroup.char[k]])
    # Add interaction test for heterogeneity
    result <- rbind(result,"")
    result[nrow(result),1] <- if (!missing(labels)) {if (length(labels) == length(attr(terms(subgroup.model),"term.labels"))) labels[k] else subgroup.char[k]} else subgroup.char[k]
    ia.pval <- anova(coxph(ia.model,data=data),coxph(main.model,data=data),test="Chisq")[2,"P(>|Chi|)"]
    result[nrow(result),ncol(result)] <- format.pval(ia.pval)
    # Add results for each subgroup level
    for (j in 1:length(factor.levels)){
      result <- rbind(result,"")
      result[nrow(result),1] <- paste("-",factor.levels[j])
      d.subgroup <- subset(data,.subgroup.var==factor.levels[j])
      result[nrow(result),2:(ncol(result)-1)] <- surv.comparison(model=base.model,data=d.subgroup,...)[3,]
    }
  }
  result
}



kmplot_t <- function(fit, rhs, main = ""){
  par(mar=c(8,9,4.1,2.1))
  tp <- c(0,30,60,90,120,180,240)
  lev <- levels(rhs)
  color <- c("#1b9e77","#d95f02","#7570b3","#e7298a")
  plot(fit, col = color, lwd = 2,
        axes = 0,  ylim=c(0,1),xlim=c(-1,280), main = main)
  at.risk <- summary(fit , times=tp, extend=T)$n.risk
  axis(1,at=seq(0,270,by=30),labels=NA)
  axis(1,at=tp)
  axis(2,at=c(0,0.25,0.5,0.75,1),pos=0)
  
  mtext(text="Survival probability",side=2,at=0.5, line=2,cex=1.2)
  mtext(text="Days since randomization",side=1,at=135, line=2,cex=1.2)
  mtext(text="No. at risk",side=1,at=-105, line=3,cex=1,adj=0)
  
  mtext(text=lev[1],side=1,at=-105,line=4,cex=0.8,adj=0)
  mtext(text=lev[2],side=1,at=-105,line=5,cex=0.8,adj=0)
  mtext(text=lev[3],side=1,at=-105,line=6,cex=0.8,adj=0)
  mtext(text=lev[4],side=1,at=-105,line=7,cex=0.8,adj=0)
  
  mtext(text=at.risk[1:7],sid=1,at=tp, line=4,cex=.8,adj=NA)
  mtext(text=at.risk[8:14],sid=1,at=tp, line=5,cex=.8,adj=NA)
  mtext(text=at.risk[15:21],sid=1,at=tp, line=6,cex=.8,adj=NA)
  mtext(text=at.risk[22:28],sid=1,at=tp, line=7,cex=.8,adj=NA)
  
  legend(x=5,y=0.25,legend=paste(lev),lty="solid",lwd=2,col=color,cex=.8, bty = "n" )
  invisible(NULL)
}

onevar.func <- function(data, tab0 = tab1, variable){
  var.level = levels(variable)
  r <- NULL
  for(i in 1:length(var.level)){
    subdat <- subset(data, variable == var.level[i])
    tab <- table(subdat$arm)
    r1 <- paste(tab, " (", formatC(tab/tab0*100,digits = 0,format = "f"),"%)", sep = "")
    r <- rbind(r,r1)
  }
  r <- cbind(t(t(var.level)), r)
  return(r)
}

my_opendex_fun <- function(data, varlist, fulldata = bl.ep){
  # header will be added at the end
  tab1 <- table(fulldata$arm)
  r0 <- paste("N = ", tab1, sep = "" )
  r0 <- c(" ",r0) 
  result <- r0
  
  for(i in 1:ncol(varlist)){
    result.i <- onevar.func(data = data, tab0 = tab1,variable = varlist[,i])
    result <- rbind(result, result.i)
  }
  rownames(result) <- NULL
  colnames(result) <- result[1,]
  result <- as.data.frame(result[-1,])
   return(result)
}

mySummary.ae <- function(ae,pt.arm,caption=NULL,filename=NA){
  require(xtable)
  # debug
  gr.lev <- levels(pt.arm$arm)
  n.1 <- sum(pt.arm$arm==gr.lev[1])
  n.2 <- sum(pt.arm$arm==gr.lev[2])
  # add dummy rows to generate numbers of AEs of any type
  ae.any <- ae; ae.any$aeterm <- 0; ae.any$aename <- "Any selected AE"
  ae <- rbind(ae,ae.any)
  # add randomized arm to AE
  ae.arm <- merge(pt.arm,ae,by="usubjid")
  
  # number of AE of each type
  # ae.arm$aename[ae.arm$aeterm==10] <- ae.arm$othaespec[ae.arm$aeterm==10]
  all.aenames <- unique(ae.arm$aename[order(ae.arm$aeterm)])
  ae.count <- data.frame("ae.name"=all.aenames,"n.pt.1"=NA,"n.ae.1"=NA,"n.pt.2"=NA,"n.ae.2"=NA,"p.val"=NA,stringsAsFactors=F)
  for (i in 1:nrow(ae.count)){
    pt.ae.1 <- ae.arm$usubjid[(ae.arm$aename==ae.count$ae.name[i])&(ae.arm$arm==gr.lev[1])]
    ae.count$n.ae.1[i] <- length(pt.ae.1)
    ae.count$n.pt.1[i] <- length(unique(pt.ae.1))
    pt.ae.2 <- ae.arm$usubjid[(ae.arm$aename==ae.count$ae.name[i])&(ae.arm$arm==gr.lev[2])]
    ae.count$n.ae.2[i] <- length(pt.ae.2)
    ae.count$n.pt.2[i] <- length(unique(pt.ae.2))
    ae.count$p.val[i] <- fisher.test(cbind(c(ae.count$n.pt.1[i],n.1-ae.count$n.pt.1[i]),c(ae.count$n.pt.2[i],n.2-ae.count$n.pt.2[i])))$p.value
  }
  ae.count.final <- cbind("Adverse event name"=ae.count$ae.name,
                          "n.pt.1"=paste(ae.count$n.pt.1," (",round(100*ae.count$n.pt.1/n.1,0),"%)",sep=""),
                          "n.ae.1"=as.character(ae.count$n.ae.1),
                          "n.pt.2"=paste(ae.count$n.pt.2," (",round(100*ae.count$n.pt.2/n.2,0),"%)",sep=""),
                          "n.ae.2"=as.character(ae.count$n.ae.2),
                          " (p value)"=as.character(round(ae.count$p.val,3)))
  ae.count.final <- rbind(c("",paste(gr.lev[1]," (n=",n.1,")",sep=""),"",paste(gr.lev[2]," (n=",n.2,")",sep=""),"","Comparison"),
                          colnames(ae.count.final),ae.count.final)                             
  if (!is.na(filename)){  # generate html table
    x <- print(xtable(ae.count.final,caption=caption),type="html",file=filename,
               caption.placement="top",include.rownames=F,include.colnames=F,html.table.attributes="border=1")
    x <- paste(x,"n.pt.1 and n.pt.2 refer to the number of patients with at least one adverse event in each study arm. n.ae.1 and n.ae.2 refer to the number of adverse events in each study arm.")
    write(x,file=filename)
  }
  ae.count.final
}

qlabel <- function(data, labels){
  vars <- colnames(data)
  if (!length(names(labels))) names(labels) <- vars
  for (var in vars)
    data[,var] <- structure(data[,var], label = labels[[var]])
  return(data)
}

`qlabel<-` <- function(x, value) qlabel(x, value)

mySummary.sort <- function(myTable, ..., marginTop = NULL, marginBottom = NULL) {
  table.sort <- if (length(marginTop) | length(marginBottom)) 
    dplyr::arrange(as.data.frame(myTable[-c(marginTop, marginBottom),]), ...)
  else
    dplyr::arrange(as.data.frame(myTable), ...)
  table.out <- rbind(myTable[marginTop,], table.sort, myTable[marginBottom,])
  rownames(table.out) <- NULL
  as.matrix(table.out)
}

mySummary.ae.sort <- function(aeTable, ..., marginTop = NULL, marginBottom = NULL){
  table.split <- aeTable[-c(marginTop, marginBottom),] 
  v.all <- as.numeric(table.split[, 3]) + as.numeric(table.split[, 5])
  table.split <- cbind(table.split, c(v.all = v.all))
  table.out <- mySummary.sort(table.split, desc(v.all), ...)[, -7]
  rbind(aeTable[marginTop,], table.out, aeTable[marginBottom,])
}
