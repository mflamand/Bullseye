
# color palette

cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# theme

library(aod)

theme_mf <- function(){
  
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color="black",size=0.25),
        axis.ticks = element_line(color="black",size=0.25),
        axis.text = element_text(color="black",size=7, margin=margin(0,0,0,0)) 
  )
  
}

# glm helper function

glm_comp_bbn<-function(df, design, lev,comp){
  
  facts_comp<-levels(df[[comp]])[[1]]
  df<-subset(df, df$cov != 0)
  df$non_mut<-df$cov-df$mut

  facts<-levels(droplevels(df[[comp]]))
  test_res<-NULL
  
  design1<-design[[2]]
  lev=paste0(design1,lev)
  
  
  if (length(facts) < 2  || !facts_comp %in% facts){ 
    p<-as.vector(rep(NA,length(lev)),mode="list")
    names(p)<-lev
   
  }
  else{
    suppressWarnings(mod<-aod::betabin(formula(paste0("cbind(mut,cov-mut)",deparse(design))), random = ~1, data=df))
    test_res<-try(vcov(mod),silent=TRUE)
  }
  
  if(is.null(test_res) || class(test_res)=="try-error"){
    p<-as.vector(rep(NA,length(lev)),mode="list")
    names(p)<-lev
    return(p)
  }
  else{
    p<-vector("list", length(lev))
    n<-vector("list", length(lev))
    LFC<-as.vector(rep(0,length(lev)),mode="list")
    names(LFC)<-levLFC
    
    okbeta <- !is.na(coef(mod, na.rm = FALSE))
    aa <- coef(mod)[okbeta]
    
    index <- which(names(aa) %in% lev)

    for(x in index){
      
      b<-coef(mod)
      V<-vcov(mod)
      L <- matrix(rep(0, length(b)), ncol = length(b))
      L[x] <- 1
      # dimnames(L) <- list(paste("L", as.character(seq(NROW(L))),sep = ""), names(b))
      H0 <- rep(0, nrow(L))
      f <- L %*% b
      
      mat <- qr.solve(L %*% V %*% t(L))
      stat <- t(f - H0) %*% mat %*% (f - H0)
      p[x] <- 1 - pchisq(stat, df = nrow(L))
      n[x] <- names(b)[x]
    }

    names(p)<- n
    
    for ( x in lev) {
      if (!x %in% names(p)){
        # print(x)
        p[[x]]<- NA   # if it was missing because of a lack of coverage, return NA as p-value
      }
    }
  }  
  return(p[ names(p) %in% lev])
}


glm_comp_bn<-function(df, design, lev,comp,link="binomial"){
  
  facts_comp<-levels(df[[comp]])[[1]]
  df<-subset(df, df$cov != 0)
  canonicalOrder <- function(term) {
    tt <- strsplit(term, ":")
    tt <- lapply(tt, sort)
    sapply(tt, paste, collapse = ":")
  }
  
  explicit1<-function (formula) {
    if (length(formula) == 1) 
      return(formula == 1)
    if (!(formula[[1]] == "+" || formula[[1]] == "*" || formula[[1]] == 
          "/" || formula[[1]] == "^" || formula[[1]] == "~")) 
      return(FALSE)
    if (length(formula) == 3) {
      (formula[[2]] == 1) || explicit1(formula[[2]]) || explicit1(formula[[3]])
    }
    else {
      (formula[[2]] == 1) || explicit1(formula[[2]])
    }
  }
  
  facts<-levels(droplevels(df[[comp]]))
  if (length(facts) < 2  || !facts_comp %in% facts){ 
    p<-as.vector(rep(1,length(lev)),mode="list")
    LFC<-as.vector(rep(0,length(lev)),mode="list")
    design<-design[[3]]
    lev=paste0(design,lev)
    levLFC<-paste0(lev,"LFC")
    names(p)<-lev
    names(LFC)<-levLFC
    lev<-c(lev,levLFC)
    p<-c(p,LFC)
    return(p[ names(p) %in% lev])
  }
  else{
    mod<-glm(formula=design,data=df,weights = cov, family=link)
    
    design<-design[[3]]
    lev=paste0(design,lev)
    levLFC<-paste0(lev,"LFC")

    p<-vector("list", length(lev))
    LFC<-as.vector(rep(0,length(lev)),mode="list")
    names(LFC)<-levLFC
    n<-vector("list", length(lev))

    design<-formula(paste0("~",design))
    
    if (inherits(design, "formula")) {
      test_intercept <- explicit1(design)
      design <- attr(terms(design), "term.labels")
    }
    else test_intercept <- FALSE
    
    okbeta <- !is.na(coef(mod, na.rm = FALSE))
    tt <- attr(terms(mod), "term.labels")
    aa <- attr(model.matrix(mod), "assign")[okbeta]
    
    index <- which(aa %in% match(canonicalOrder(design), canonicalOrder(tt)))
    
    if (test_intercept) {
      if (attr(terms(model), "intercept")) 
        index <- unique(c(1, index))
      else stop("model does not have an intercept")
    }
    
    # return(index)
    for(x in index){
      beta <- coef(mod)[x]
      name_fact<-(names(beta))
      V<-NULL
      # V <- vcov(mod)[x, x]
      V<-try(vcov(mod)[x, x],silent=TRUE)
      if(is.null(V) || is.na(V) || class(V)=="try-error" || V == 0){
        p[x]<-1
      }
      else{
        chisq <- beta %*% solve(V) %*% beta
        pval = pchisq(chisq, length(2), lower.tail = FALSE)
        p[x]<-pval[1,1] # get p-value from result
      }
      n[x]<-name_fact
      LFC[x]<-beta
    }
    names(p)<- n
    # names(LFC)<- paste0(n,"LFC")
    
    for ( x in lev) {
      if (!x %in% names(p)){
        p[[x]]<- NA   # if it was missing because of a lack of coverage, return NA as p-value
        LFC[[paste0(x,"LFC")]]<- 0
      }
    }
  } 
  p<-c(p,LFC)
  lev<-c(lev,levLFC)
  return(p[ names(p) %in% lev])
  
}


#function to calculate glm and p value

Bullseye<-function(data_cov,data_mut, colData=NULL, filter.col=NULL, filter.with=NULL, min=0.05, filter.sites=NULL,dataset=NULL,design=NULL, min.cov=10, return.df.only=FALSE, link='betabin', numCores=1){

  #filter.sites : needs to be "all" or "called", if all: will only keep sites with at least min edit in one or the other samples. If "called": will keep samples with at least min in both conditions. 
  if(link != 'betabin' && link != 'bin' && link != 'quasibinomial') stop('Please use link="betabin", link="bin" or link="quasibinomial')
  
  if(is.null(colData)) stop("please provide colData with sample information for filtering")
  if(is.null(design)) stop("please provide design formula with a factor found in colData (eg. ~Genotype)")
  
  if(!is.null(filter.col)){
   colData<- colData[colData[[filter.col]]==filter.with,] #base version 
  }
  
  if(!is.factor(colData[[ design[[2]] ]])){
    colData[[ design[[2]] ]]<- as.factor(colData[[ design[[2]] ]])
  }
  
  colData[[ design[[2]] ]]<- droplevels(colData[[ design[[2]] ]])
  
  design_levels<-levels(colData[[ design[[2]] ]])

  #only keep matrix columns from columns in colData 
 
  if(! "n" %in% colnames(data_cov))
  {
    data_cov$n = 1
    data_mut$n = 1
  }

  data_cov<-data_cov[,c("n",rownames(colData))]
  data_mut<-data_mut[,c("n",rownames(colData))]
 
  # We can filter based on a list of sites found in .dataset. Else we can filter and keep a subset of samples from colData.
  if(!is.null(dataset)){
    dataset$site<- paste(dataset$chr,dataset$start,dataset$end,dataset$strand,sep="_")
    data_cov<-data_cov[rownames(data_cov) %in% dataset$site,]
    data_mut<-data_mut[rownames(data_mut) %in% dataset$site,]
  }
  
  n_site<-data_cov$n

  # filter based on min coverage, replace with 0 to exclude from analysis
  data_cov[data_cov < min.cov*data_cov$n] <- 0
  data_mut[data_cov < min.cov*data_mut$n] <- 0

  data_cov$n<-n_site
  data_mut$n<-n_site
  # calculate editing rate in each sample
  data_edit<-(data_mut/data_cov)*n_site
  data_edit[is.na(data_edit)]<-NA
  
  # first we get the average editing in called sites
  data_edit_perc<-data_edit
  data_edit_perc[data_edit_perc< min]<-NA
  data_edit_perc[is.na(data_edit_perc)]<-NA  # NaN to NA if any

  data_cov <- data_cov %>% rownames_to_column("site") # change to tibble style for pivoting and nesting of dataframe
  data_mut <- data_mut %>% rownames_to_column("site")
  data_edit <- data_edit %>% rownames_to_column("site")
  data_edit_perc <- data_edit_perc %>% rownames_to_column("site")

  L_cov<- data_cov %>%
    pivot_longer(names_to = "sample",cols = !c(site,n) ,values_to = "cov") %>%
    left_join(., colData %>% rownames_to_column("sample"), by="sample")

  L_mut<-data_mut %>%
    pivot_longer(names_to = "sample",cols = !c(site,n),values_to = "mut")
  
  data_long <- suppressMessages(dplyr::left_join(L_cov,L_mut))
  n_data<-data_long%>% group_by(site) %>% nest()
  
  if (link =="bin"){
    results<-parallel::mclapply(X = n_data$data, FUN = glm_comp_bn,design=formula(paste0("mut/cov",deparse(design))), lev=design_levels,comp=design[[2]], mc.cores = numCores)
  }
  else if (link =="quasibinomial"){
    results<-parallel::mclapply(X = n_data$data, FUN = glm_comp_bn,design=formula(paste0("mut/cov",deparse(design))), lev=design_levels,comp=design[[2]],link="quasibinomial", mc.cores = numCores)
  }
  else if (link =="betabin"){
    results<-parallel::lapply(X = n_data$data, FUN = glm_comp_bbn,design=design, lev=design_levels,comp=design[[2]], mc.cores = numCores)
  }
  # return(results)
  results<-do.call(rbind.data.frame, results)
  # return(results)
  site<- n_data$site
  results<-cbind(site,results)
 
  # first level is always first and is fixed.
  first.level<-design_levels[[1]]
  first.level.regex<-paste0(".*",first.level,".*")
  
  data_edit_A<-data.frame("site"='')
  
  i<-length(design_levels)
  while( i >= 2){
    
    second.level<-design_levels[[i]]
    second.level.regex<-paste0(".*",second.level,".*")

    data_edit_A<-suppressMessages(dplyr::mutate(.data=data_cov,
                                                !!paste0(first.level,"_sum_cov") := rowSums(select(data_cov, matches( first.level.regex,perl=TRUE)),na.rm=TRUE)) %>% 
                                    dplyr::mutate(!!paste0(first.level,"_mut") := rowSums(select(data_mut, matches(first.level.regex,perl=TRUE)),na.rm=TRUE)) %>% 
                                    mutate(!!first.level := .data[[paste0(first.level,"_mut")]]/.data[[paste0(first.level,"_sum_cov")]]) %>% 
                                    mutate(!!paste0(first.level,"_min_cov") := rowSums(!is.na(select(.data=data_edit, matches(first.level.regex,perl=TRUE)))),
                                           !!paste0(first.level,"_called_sites") := rowSums((select(.data=data_edit, matches(first.level.regex,perl=TRUE)))   >= min,na.rm=TRUE )) %>% 
                                    mutate(!!paste0(first.level,"_average_called") := rowMeans(select(.data=data_edit_perc, matches(first.level.regex,perl=TRUE)),na.rm = TRUE)) %>% 
                                    dplyr::mutate(!!paste0(second.level,"_sum_cov") := rowSums(select(.data=data_cov, matches( second.level.regex,perl=TRUE)),na.rm = TRUE)) %>%
                                    dplyr::mutate(!!paste0(second.level,"_mut") := rowSums(select(.data=data_mut, matches(second.level.regex,perl=TRUE)),na.rm = TRUE)) %>% 
                                    mutate(!!second.level := .data[[paste0(second.level,"_mut")]]/.data[[paste0(second.level,"_sum_cov")]],
                                           !!paste0(second.level,"_log2FoldChange") := log2(.data[[second.level]]/.data[[first.level]])) %>% 
                                    mutate(!!paste0(second.level,"_min_cov") := rowSums(!is.na(select(.data=data_edit, matches(second.level.regex,perl=TRUE)))),
                                           !!paste0(second.level,"_called_sites") := rowSums((select(.data=data_edit, matches(second.level.regex,perl=TRUE))) >= min,na.rm=TRUE)) %>% 
                                    mutate(!!paste0(second.level,"_average_called") := rowMeans(select(.data=data_edit_perc, matches(second.level.regex,perl=TRUE)),na.rm = TRUE)) %>% 
                                    left_join(.,data_edit_A) %>% left_join(.,results %>% select(site, !!paste0(second.level,"_pvalue"):=matches( paste0(".*",second.level,"$")), 
                                                                                                !!paste0(second.level,"_LFC"):=matches(paste0(".*",second.level,"LFC") )) ) )
    
    i<-i-1
  }
  
  data_edit_A<- data_edit_A %>% mutate(baseEdit = rowSums(select(., matches("average_called")),na.rm = TRUE)*100*n/rowSums(select(., matches("called_sites")),na.rm = TRUE)) %>% mutate(baseEdit=replace_na(baseEdit,0))
  data_edit_A<-data_edit_A %>% dplyr::select(-rownames(colData))
  
  
  if (!is.null(filter.sites)){
    R.utils::printf("Filtering sites with less than %.1f %% editing", min*100)
    if (filter.sites == "called"){
      data_edit_A<-data_edit_A %>%
        dplyr::filter(baseEdit >= min)
      # dplyr::filter(.data[[paste0(first.level,"_average_called")]] >= min & .data[[paste0(second.level,"_average_called")]] >= min)
    }
    else if (filter.sites == "all"){
      data_edit_A<-data_edit_A %>%
        dplyr::filter(.data[[first.level]] >= min | .data[[second.level]] >= min)
    }
  }
  
  return(data_edit_A)
  
}



