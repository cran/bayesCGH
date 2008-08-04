setClass("ClassifySegList", representation("SegList"),
         prototype = prototype(list(state = NULL, rpred = NULL, prob = NULL, M.predicted = NULL, dispersion = NULL, variance = NULL, M.observed = NULL, genes = NULL, probes.class = "matrix", regions = "list")),
         contains = "SegList")

setMethod("dim", signature = "ClassifySegList", function(x) {callNextMethod(x)})
setMethod("length", signature = "ClassifySegList", function(x) {callNextMethod(x)})
setMethod("dimnames", signature = "ClassifySegList", function(x) {callNextMethod(x)})

setMethod("show", signature = "ClassifySegList", function(object) {
	cat("An object of class \"",class(object),"\"\n",sep="")
	for (what in setdiff(names(object), c("regions", "probes.class"))) {
		y <- object[[what]]
		cat("$",what,"\n",sep="")
		printHead(y)
		cat("\n")
	}
	for (what in setdiff(slotNames(object), c("regions", "probes.class", ".Data"))) {
		y <- slot(object,what)
		if(length(y) > 0) {
			cat("@",what,"\n",sep="")
			printHead(y)
			cat("\n")
		}
	}
})


