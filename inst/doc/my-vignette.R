## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=FALSE, echo=T, results='hide',message=F, warning=F-----------------
#  install_github("LuisBenavides/SNSchart")

## ---- echo=TRUE, eval=FALSE, results='hide'-----------------------------------
#  library("SNSchart")

## ---- echo=TRUE, eval=FALSE, results='hide'-----------------------------------
#  X = SNSchart::example82$X
#  X.id = SNSchart::example82$X.id

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(head(SNSchart::example82, 10))

## ----echo=T,eval=F, results='hide', message=F, warning=F----------------------
#  s = SNSchart::SNS(X=X,X.id=X.id)

## ----echo=T,eval=F, results='hide', message=F, warning=F----------------------
#  plot(s)

## ----echo=F,eval=T, results='hide', message=F, warning=F,fig.width = 9, fig.height = 6----
X = SNSchart::example82$X #get the dataset into a data frame
X.id = SNSchart::example82$X.id
s = SNSchart::SNS(X=X,X.id=X.id)
plot(s)

## ---- echo=TRUE, eval=FALSE, results='hide'-----------------------------------
#  X = SNSchart::example84$X
#  X.id = SNSchart::example84$X.id

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(head(SNSchart::example84, 10))

## ----echo=T,eval=F------------------------------------------------------------
#  s = SNSchart::SNS(X=X,X.id=X.id, chart="CUSUM", chart.par=c(0.5, 4.389, 3))

## ----echo=T,eval=F------------------------------------------------------------
#  plot(s)

## ----echo=F,eval=T, results='hide', message=F, warning=F,fig.width = 9, fig.height = 6----
X = SNSchart::example84$X #get the dataset into a data frame
X.id = SNSchart::example84$X.id
s = SNSchart::SNS(X=X,X.id=X.id, chart="CUSUM", chart.par=c(0.5, 4.389, 3)) 
plot(s)

## ---- echo=TRUE, eval=FALSE, results='hide'-----------------------------------
#  X = SNSchart::example84$X
#  X.id = SNSchart::example84$X.id

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(head(SNSchart::example84, 10))

## ----echo=T,eval=F------------------------------------------------------------
#  s = SNSchart::SNS(X=X,X.id=X.id, chart="EWMA", chart.par=c(0.01, 2.0171))

## ----echo=T,eval=F------------------------------------------------------------
#  plot(s)

## ----echo=F,eval=T, results='hide', message=F, warning=F,fig.width = 9, fig.height = 6----
X = SNSchart::example84$X #get the dataset into a data frame
X.id = SNSchart::example84$X.id
s = SNSchart::SNS(X=X,X.id=X.id, chart="EWMA", chart.par=c(0.01, 2.0171))
plot(s)

## ---- echo=TRUE, eval=FALSE, results='hide'-----------------------------------
#  X = SNSchart::example87$X
#  X.id = SNSchart::example87$X.id
#  Y = SNSchart::example87$Y

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(head(SNSchart::example87, 10))

## ----echo=T,eval=F------------------------------------------------------------
#  s = SNSchart::SNS(X=X,X.id=X.id, Y=Y, chart="EWMA", chart.par=c(0.01, 2.0171))

## ----echo=T,eval=F------------------------------------------------------------
#  plot(s)

## ----echo=F,eval=T, results='hide', message=F, warning=F,fig.width = 9, fig.height = 6----
X = SNSchart::example87$X #get the dataset into a data frame
X.id = SNSchart::example87$X.id
Y = SNSchart::example87$Y
s = SNSchart::SNS(X=X,X.id=X.id, Y=Y, chart="EWMA", chart.par=c(0.01, 2.0171))
plot(s)

## ---- echo=TRUE, eval=FALSE, results='hide'-----------------------------------
#  X = example49$X2
#  X.id = example49$X.id
#  Y = example49$Y2

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(head(SNSchart::example49, 10))

## ----echo=T,eval=F------------------------------------------------------------
#  s = SNSchart::SNS(X=X,X.id=X.id, Y=Y, chart="Shewhart", scoring="Z-SQ",isFixed = TRUE)

## ----echo=T,eval=F------------------------------------------------------------
#  plot(s)

## ----echo=F,eval=T, results='hide', message=F, warning=F,fig.width = 9, fig.height = 6----
X = SNSchart::example49$X2
X.id = SNSchart::example49$X.id
Y = SNSchart::example49$Y2
s = SNSchart::SNS(X=X,X.id=X.id, Y=Y, chart="Shewhart", scoring="Z-SQ",isFixed = TRUE)
plot(s)

## ---- echo=TRUE, eval=FALSE, results='hide'-----------------------------------
#  X = SNSchart::example91[,1:2]
#  X.id = SNSchart::example91$X.id

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(head(SNSchart::example91, 10))

## ----echo=T,eval=F------------------------------------------------------------
#  msns = SNSchart::MSNS(X, X.id)

## ----echo=T,eval=F------------------------------------------------------------
#  plot(msns)

## ----echo=F,eval=T, results='hide', message=F, warning=F,fig.width = 9, fig.height = 6----
X = SNSchart::example91[,1:2] #get the dataset into a data frame
X.id = SNSchart::example91$X.id
msns = SNSchart::MSNS(X, X.id)
plot(msns)

## ---- echo=TRUE, eval=FALSE, results='hide'-----------------------------------
#  X = SNSchart::example93[,1:2]
#  X.id = SNSchart::example93$X.id

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(head(SNSchart::example93, 10))

## ----echo=T,eval=F------------------------------------------------------------
#  msns = SNSchart::MSNS(X, X.id, null.dist = "F")

## ----echo=T,eval=F------------------------------------------------------------
#  plot(msns)

## ----echo=F,eval=T, results='hide', message=F, warning=F,fig.width = 9, fig.height = 6----
X = SNSchart::example93[,1:2] #get the dataset into a data frame
X.id = SNSchart::example93$X.id
msns = SNSchart::MSNS(X, X.id, null.dist = "F")
plot(msns)

