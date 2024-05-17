#function to be used for forward selection ordinal
# pw.auc <- function(probs,true){
#     auc0=AUC(probs[true$dhy_ord==0,1],probs[true$dhy_ord!=0,1])
#     auc2= AUC(probs[true$dhy_ord==2,3],probs[true$dhy_ord!=2,3])
#     pw.auc <- (auc0 + auc2)/2
#     return(pw.auc)
# }
library(matlab)
rcs.fml.fnc <- function(splines=T,variable,data,cont.vars){
  if(variable %in% cont.vars && splines==T){
    quantile.knots = paste('c(',
                           paste(quantile(data[,variable], c(0.1, 0.5, 0.9), 
                                          include.lowest = TRUE,na.rm=T),collapse=','),')')
    rcs.fml <- paste('rcs(',variable,',',quantile.knots,')')
  }else{
    rcs.fml <- variable
  }
  return(rcs.fml)
}
interaction.pairs <- function(j,predictors.in){
  indicies <- expand.grid(1:j, 1:j)[expand.grid(1:j, 1:j)[,1]!=expand.grid(1:j, 1:j)[,2],]
  interact.predictors <- as.matrix(data.frame(t(apply(indicies,1,function(x){predictors.in[sort(x)]}))) %>% distinct())
  return((interact.predictors))
}
pdispecific=function(y,probs,class,ysep,c,k){
    other=1:k;
    other=other[-class]
    matrix=matrix(0,prod(c),k)
    matrix[,1]=repmat(t(probs[ysep[class]$y,class]),prod(c[-class]),1)
    
    for (i in 1:(length(other)-1)){
        matrix[,1+i]=repmat(
        t(sort( repmat(t(probs[ysep[other[i]]$y,class]),c[class]*max(1,prod(c[other[1:i-1]])),1) )),prod(c[other[(i+1):length(other)]]),1)
    }
    i=length(other)
    matrix[,i+1]=t(sort( repmat(t(probs[ysep[other[i]]$y,class]),c[class]*max(1,prod(c[other[1:i-1]])),1) ))
    # there are no ties when the last index of the max within each row equals 1 -> count the number of rows for which this occurs
    pdispecific=sum(max.col(matrix,ties.method="last")==1) # number of cases for which probs is max for case of class i
    # handle ties: the number of sets with ties is defined as the number of row in the matrix for which the last index of the max value within that row differs from 1 AND the first index equals 1 (otherwise max chance not for case of the correct class)
    temp =(matrix[max.col(matrix,ties.method="last")!=1&max.col(matrix,ties.method="first")==1,]- matrix[max.col(matrix,ties.method="last")!=1&max.col(matrix,ties.method="first")==1,1]==0)*1
    tievalue=sum(1/sum(t(temp)))
    pdispecific=(pdispecific+tievalue)/nrow(matrix)
}


pdivalue=function(y,probs){
    library(matlab)
    pdi=0
    w=table(y)
    c=as.numeric(w) # number of cases per category
    k=length(c) # number of categories
    # get for each category the location of its cases in the dataset
    ysep=c()
    for (i in 1:k) {
        ysep=c(ysep,list(y=which(y==i)))
    }
    # calculate class-specific pdi_i for each class i
    pdi_i=rep(0,k)
    for (i in 1:k){ pdi_i[i]=pdispecific(y,probs,i,ysep,c,k) }
    pdivalue=c(sum(pdi_i)/k,pdi_i)
}


auc_bin <- function(probs,true.label){
    auc01=AUC(probs,true.label)
    return(auc01)
}

ManualCalculateProbs <- function(model,formula,train){
    mm <- model.matrix(formula,data= train) #expand factor to a set of dummy variables and expand interactions 
    #some coefficients dropped due to rank deficiency
    colnms.in <- which(colnames(mm) %in% names(model$coefficients))# r
    mm <- mm[,colnms.in]
    preds <- matrix(ncol=1,nrow=nrow(mm))
    sub.idx = which(!is.na(coef(model)))
    preds <- mm[,sub.idx]%*% coef(model)[sub.idx] 
    preds <- 1/(1+exp(-preds))
    return(preds)
}
ManualCalculatelglkBin = function(label.true,probs){
  label.true = match(label.true, sort(unique(label.true)))
  need.info = cbind(label.true, 1-probs,probs)
  lglk = sum( apply(need.info,1,function(x){
        log(x[-1][x[1]])}))
         return(lglk)
}

#function to implement parallel computing: it iterates through each pair of all candidate interactions in order to find
#the best pair in the next step ; returns the evalution index: AIC, auc, or misclassification rate.
iterateLeftPairsBin <- function(response,pairs,interaction.formula,train,eval.method){
    pairs  <- paste(pairs, collapse="*")
    interaction.formula.j <- update(interaction.formula, as.formula(paste('.~ .+',pairs)))
    err_catch <- safely(function(formula) {glm(formula, family="binomial", data = train)})
    if(!is.null(err_catch(interaction.formula.j)$error)){
        index.interact.j <- NA
    } else{
        model <- glm(interaction.formula.j, family="binomial", data=train)
        probs <- ManualCalculateProbs(model,interaction.formula.j,train)
        #the following ifelse takes account for rank deficient model, coefs might be dropped off
        if(eval.method=='lglk'){
          index.interact.j <- -logLik(model)
        }else if (eval.method=='aic'){
            index.interact.j <-AIC(model)
        } else if(eval.method=='mindex'){
            mindex <- auc_bin(probs, train[,response])
            index.interact.j <-  -mindex#make negative so best performance has the min index
        }else if(eval.method=='pdi'){
            #print(interaction.formula.j)
            pdivalues <- pdivalue(as.numeric(train[,response]),probs)
            index.interact.j <- -pdivalues[1]
        }else{
            pdivalues <- pdivalue(as.numeric(train[,response]),probs)
            index.interact.j <- -sum(c(table(train[,response])/nrow(train)*pdivalues[-1]))
        }
    }
    #print(interaction.formula.j)
    #print(index.interact.j)
    return(index.interact.j)
}


#use results from iterateLeftPairsBin to find the best interaction pair
#for each interaction size of a main size, do this:
find.interactionsBin<- function(predictor.format,response,train,test,eval.method,ncore=2){
  cl2 <- makeCluster(ncore, type="FORK")
  registerDoParallel(cl2)
  orc.interact.i <- c()
    if(sum(is.na(predictor.format))==0){
        interaction.formula <- as.formula(paste(response, paste(predictor.format, collapse=" + "), sep=" ~ "))
        i=length(predictor.format)
        combo.n <- choose(i, 2)
        combos.out <-data.frame(interaction.pairs(i,predictor.format),stringsAsFactors=FALSE)
        combos.in <-c()
        for ( h in 1:combo.n){
            index.interact <-c()
         
            index.interact <- foreach (j=1:nrow(combos.out),
                                       .combine = c,
                                       .packages = c("doParallel")) %dopar% {
                                         iterateLeftPairsBin(response,combos.out[j,],interaction.formula,train,eval.method)
                                       } #change %dopar% to %do% if not run parallel
            
            if(length(which.min(index.interact))==0){
                print('break')
                orc.interact.i[h:combo.n] <-rep(NA,combo.n-h+1)
                break # break out of h loop
            } else{
              interaction.formula.upd <- update(interaction.formula,
                                                  as.formula(paste('.~ .+',
                                                                   paste(combos.out[which.min(index.interact),], collapse="*"))))
                # #interaction.format.i[[h]]<- paste(combos.out[which.min(index.interact),], collapse="*")
                combos.out <- combos.out[-which.min(index.interact),] #update combos.out
                interaction.model.upd <- glm(interaction.formula.upd,family="binomial",train)
                probs.test <- ManualCalculateProbs(interaction.model.upd,interaction.formula.upd,test)
                #probs.test <- predict(interaction.model.upd,test,type='probs')
                lglk.test = ManualCalculatelglkBin(test[,response],probs.test)
                # preds.ord <- apply(probs,1,which.max)-1
                orc.interact.i[h] <- -lglk.test
            }
        }
        stopCluster(cl2)
    } else{
        orc.interact.i <- rep(NA,choose(length(predictor.format),2) )
    }
    return( orc = orc.interact.i)
}

#function to find the best main size and best interaction size
FindModelSizeBin <- function(k,response,names,splines=F,include.interaction,
train.data,eval.method){
    err.msg <- err.model <- orc.interact <- model.formulas <-predictor.format.list<- list()
    predictors.in <-orc.main <- orc.interact.i <- c()
    p  <- length(names)
    train <- train.data[train.data$folds!=k,]
    test <- train.data[train.data$folds ==k,]
    for (i in 1:p){
        index.aic <- mindex <- pdi <- pdi.wt <- lglk <- c()
        predictors.out <- names[!names %in% predictors.in] #predicotrs not selected in model yet
        #main effect#
        for (j in 1:length(predictors.out)){
            predictors <- c(predictors.in, predictors.out[j])
            predictor.format <- sapply(predictors,function(x){rcs.fml.fnc(splines,x,train,cont.vars)})
            model.formula <- as.formula(paste(response, paste(predictor.format, collapse=" + "), sep=" ~ "))
            #skip and catch errors
            err_catch <- safely(function(formula) {glm(formula,family="binomial",data = train)})
            if(!is.null(err_catch(model.formula)$error)){
                err.model <- list(err.model,model.formula)
                err.msg <- list(err.msg,err_catch(model.formula)$error )
                index[j] <- NA
            } else {
                model <- glm(model.formula, family="binomial", data=train)
                mm <- model.matrix(model.formula,train[,c(response,names)])
                #the following ifelse takes account for rank deficient model, coefs might be dropped off
                if(ncol(mm) != length(coef(model))) {#different if for oridnal models
                    index.aic[j] <- mindex[j] <- pdi[j] <- pdi.wt[j] <- lglk[j] <- NA
                }else{
                    lglk[j] <- -logLik(model)
                    index.aic[j] <- AIC(model)
                    probs <- predict(model, train[,-1],type='response')
                    mindex[j] <- -auc_bin(probs, train[,response])
                    if(eval.method %in% c('pdi','pdi.wt')){
                        pdivalues <- pdivalue(as.numeric(train[,response],probs))
                        pdi.wt[j] <- -sum(c(table(train[,response])/nrow(train)*pdivalues[-1]))
                        pdi[j] <- -pdivalues[1]
                    }
                }
            }
        }
        if (eval.method=='lglk'){
          index <- lglk
        } else if(eval.method=='aic'){
            index <-index.aic
        }else if(eval.method=='mindex'){
            index <-mindex
        } else if(eval.method=='pdi'){
            index <-pdi
        }else{
            index <-pdi.wt
        }
        ##if no more than two levels of a factor (not eough data, break iter. of i)
        if (length(predictors.out[which.min(index)])==0){
            orc.main[i:p] <- rep(NA,(p-i+1))
            for(b in i:p){#orc.interact[b] <- rep(NA,choose(b,2));
                predictor.format.list[[b]]  <- rep(NA,b)
            }
            break
        } else {
            predictors.in <-c(predictors.in,
            predictors.out[which.min(index)]) #min or max depends on the criteria
            predictor.format <- sapply(predictors.in,function(x){rcs.fml.fnc(splines,x,train,cont.vars)})
            model.formula.upd <- as.formula(paste(response, paste(predictor.format, collapse=" + "), sep=" ~ "))
            model.upd = glm(model.formula.upd,family = 'binomial',train)
            probs.test <- predict(model.upd,test,type='response')
            lglk.test = ManualCalculatelglkBin(test[,response],probs.test)
            orc.main[i]<- -lglk.test
            predictor.format.list[[i]] <- predictor.format
        }
    }
    
    orc.interact <- list()
    if(include.interaction ==T){
        orc.interact[[1]] <- NA
        for (i in 2:(length(names))){
            pl = predictor.format.list[[i]]
            orc.interact[[i]] <- find.interactionsBin(pl,response,train,test,eval.method)
        }
    }
    return(list(orc.main = orc.main, orc.interact = orc.interact))
}
####
#this function integrates all the functions, the steps are
#1.get the best model size (main, interaction) with k.forward cross validation;
#2. find the best predictors and interaction pairs with the just found model size on the entire training data
#3. test the model obtained from step 2 on testing data
#funtion returns best model from step 2, pairwise AUC, predicted outcomes of testing data
forward.selection.modelBinary <- function(K.forward,response,names,splines,include.interaction, 
                                           train.data,test.data,eval.method='orc'){
    predictors.in <- c()
    p <- length(names)
    #orc <- matrix(nrow = K.forward, ncol = p)
    err.msg <- err.model <- orc.interact <- orc.interact <- model.formulas <- list()
    orc.all.main <- orc.all.inta <-list()
    for (k in 1:K.forward) {

        info = FindModelSizeBin(k=k,response,
                                names=names,
                                splines=splines,
        include.interaction=include.interaction,
        train.data=train.data,eval.method)
        orc.all.main[[k]] <-info$orc.main
        orc.all.inta[[k]] <-info$orc.interact
    }
    #savenames <- paste0('ncvage',which.age,m,method)
    #save(orc.all.main,file = paste0(savenames,"ordmain.RData"))
    #save(orc.all.inta,file='/Users/xuenanwang/Desktop/Xuenan Wang thesis/Model Building/orc_all_inta.RData')
    orc.int.avg <-list();orc.int.avg[[1]] <-NA
    if(include.interaction==T){
      for (main.size in 2:p){
        amatrix <- sapply(orc.all.inta, "[[",main.size)
        #a matrix for each main sizecolmn: each fold; row: main orc +interaction orc's
        if( is.matrix(amatrix)){
          orc.int.avg[[main.size]] <- apply(amatrix, 1,
                                            function(x)mean(x,na.rm=TRUE))
        }else {
          orc.int.avg[[main.size]] <- mean(amatrix)#sapply(amatrix, function(x)mean(x,na.rm=TRUE))
        }
        rm(amatrix)
      }
      main.size.best  <- which.min(lapply(orc.int.avg,function(x) max(x)))
      #print(orc.int.avg[[main.size.best]])
      inter.size.best<- which.min(orc.int.avg[[main.size.best]])
      
    }else{
      inter.size.best= 0
      amatrix = matrix(unlist(orc.all.main), ncol = K.forward)
      main.size.best = which.min(apply(amatrix, 1,
                                       function(x)mean(x,na.rm=TRUE)))
    }
    
    
    predictors.in <- c()
    orc.main <- c()
    #main effect
    for (i in 0:(main.size.best-1)){
        index <- c()
        predictors.out <- names[!names %in% predictors.in]
        #main effect#
        index.aic <- mindex <- pdi <- pdi.wt <- lglk <- c()
        for (j in 1:length(predictors.out)){
            predictors <- c(predictors.in, predictors.out[j])
            predictor.format <- sapply(predictors,function(x){rcs.fml.fnc(splines,x,train.data,cont.vars)})
            model.formula <- as.formula(paste(response, paste(predictor.format, collapse=" + "), sep=" ~ "))
            #skip and catch errors
            err_catch <- safely(function(formula) {glm(formula, family="binomial", data=train.data)})
            if(!is.null(err_catch(model.formula)$error)){
                err.msg <- list(err.msg,err_catch(model.formula)$error )
                index[j] <- NA
            } else {
                model <- glm(model.formula, family="binomial", data=train.data)
                mm <- model.matrix(model.formula,train.data[,c(response,names)])
                if(ncol(mm) != length(coef(model))) {#more accurately ,use >
                    index.aic[j] <- mindex[j] <- pdi[j] <-pdi.wt[i]<- lglk[j] <- NA
                }else{
                    lglk[j] <- -as.numeric(logLik(model))
                    index.aic[j] <- AIC(model)
                    probs <- predict(model, train.data[,-1],type='response')
                    mindex[j] <- -auc_bin(probs,train.data[,response])
                    if(eval.method %in% c('pdi','pdi.wt')){
                        pdivalues <- pdivalue(as.numeric(train.data[,response]),probs)
                        pdi.wt[j] <- -sum(c(table(train.data[,response])/nrow(train.data)*pdivalues[-1]))
                        pdi[j] <- -pdivalues[1]
                    }
                }
            }
        }
        if(eval.method=='lglk'){
          index <- lglk
        } else if(eval.method=='aic'){
            index <-index.aic
        }else if(eval.method=='mindex'){
            index <-mindex
        } else if(eval.method=='pdi'){
            index <-pdi
        }else{
            index <-pdi.wt
        }
        predictors.in <-c(predictors.in,
        predictors.out[which.min(index)]) #min or max depends on the criteria
        orc.main[i] <-  index[which.min(index)]
    }
    predictor.format = sapply(predictors.in,function(x){rcs.fml.fnc(splines,x,train.data,cont.vars)})
    model.formula.upd <- as.formula(paste(response, paste(predictor.format, collapse=" + "), sep=" ~ "))
    model.upd <- glm(model.formula.upd, family="binomial", data=train.data)
    ###interaction####
    interaction.formula <- model.formula.upd
    orc.interact  <- c()
    if(include.interaction==F){
        orc.interact <- NA
        interaction.model <- model.upd
    } else{
        combo.n <- choose(main.size.best, 2)
        combos.out <-data.frame(interaction.pairs(main.size.best,predictor.format),stringsAsFactors=FALSE)
        for ( h in 1:inter.size.best){
            index.interact <- index.aic <- mindex <- pdi <- pdi.wt<- lglk <- c()
            #find the best pair
            for ( j in 1:nrow(combos.out)){
                pairs  <- paste(combos.out[j,], collapse="*")
                #print('before formula update')
                interact.formula <- update(interaction.formula, as.formula(paste('.~ .+',pairs))) # interact formula  =interact.formula + one pair
                #print('pass here')
                err_catch <- safely(function(formula) {glm(formula,family="binomial",data = train.data)})
                if(!is.null(err_catch(interact.formula)$error)){
                    index.interact[j] <- NA
                } else{
                    #print(!is.null(err_catch(interact.formula)$error))
                    model <- glm(interact.formula, family="binomial", data=train.data)
                    #print(coef(model))
                    probs <- ManualCalculateProbs(model,interact.formula,train.data)
                    #preds.ord <- apply(probs,1,which.max)-1
                    index.aic[j] <- AIC(model)
                    lglk[j] <- -as.numeric(logLik(model))
                    #mindex[j] <- -auc_bin(probs, train.data[,response])
                    if(eval.method %in% c('pdi','pdi.wt')){
                        pdivalues <- pdivalue(as.numeric(train.data[,response]),probs)
                        pdi.wt[j] <- -sum(c(table(train.data[,response])/nrow(train.data)*pdivalues[-1]))
                        pdi[j] <- -pdivalues[1]
                    }
                }
            }
            if(eval.method=='lglk'){
               index.interact <- lglk
            } else if (eval.method=='aic'){
                index.interact <-index.aic
            }else if(eval.method=='mindex'){
                index.interact <-mindex
            } else if(eval.method=='pdi'){
                index.interact <-pdi
            }else{
                index.interact <-pdi.wt
            }
            
            if(length(which.min(index.interact))==0){
                orc.interact[h] <-NA
                break # break out of h loop
            }
            else{
                interaction.formula <- update(interaction.formula,
                as.formula(paste('.~ .+',
                paste(combos.out[which.min(index.interact),], collapse="*"))))
                combos.out <- combos.out[-which.min(index.interact),]
                interaction.model <- glm(interaction.formula,data=train.data,family="binomial")
                orc.interact[h] <- index.interact[which.min(index.interact)]
                #print(paste('h',h))
            }
        }
    }
    
    forward.fit <- interaction.model
    forward.pred.prob <- ManualCalculateProbs(forward.fit,interaction.formula,test.data)
    #forward.orc <- mean(forward.pw.auc)
    return(list( model=forward.fit,
    pred.prob = forward.pred.prob,
    size =c(main.size.best,inter.size.best),
    stepwise.idx.main = orc.main,
    stepwise.idx.int = orc.interact))
}
