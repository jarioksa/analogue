#' Gower Distance for Mixed Data with Tied Ranks
#'
#' Proof-of-the-concept function for implementing Gower distance
#' optionally with Podani's amendment to handling ordered factors as
#' metric data with ties.
#'
#' The function accepts a data frame of potentially different data
#' types, and used in mixed metrics of Gower distance. The data types
#' for each variable are interpreted from the data frame using
#' strategies given in \code{ordered} and \code{continuous} or as
#' given for each variable in argument \code{vtypes}.
#'
#' Factors are always interpreted as nominal variables (\code{vtypes =
#' "N"}) and analysed as Matching coefficient.
#'
#' With \code{orderred = "internal"} factor levels are replaced with
#' their internal values from 1 to number of levels (\code{vtype =
#' "O"}). Strategy \code{"metric"} (\code{vtypes = "OR"}) replaces
#' them with their averaged ranks from 1 to number of observations
#' (but due to ties usually excluding these extremes). Strategy
#' \code{ordered = "podani"} (\code{vtypes = "OP"}) does the same but
#' takes into account the ties in the analysis following Podani
#' (1999). Finally, strategy \code{"nominal"} analyses them as
#' unordered factors.
#'
#' Continuous variables can be analysed as such with strategy
#' \code{continuous = "metric"} (\code{vtypes = "C"}) using Manhattan
#' distance scaled to maximum = 1. With \code{continuous = "presence"}
#' (\code{vtypes = "P"}), double-zeros are ignored for non-negative
#' variables (but variables with negative values are treated as
#' \code{"C"}), and this gives an index that is similar to a
#' quantitative variant of Jaccard dissimilarity index. Finally, with
#' strategy \code{"binary"} double-zeros are ignored in binary
#' variables but all other numeric variables are treated as
#' \code{"C"}. This gives Jaccard index to binary variables. This case
#' is known as \dQuote{asymmetric binary} in some other
#' implementations. There is no explicit \dQuote{symmetric binary}
#' because this is similar as treating variable as \code{"C"}.

#' @return Function returns a numeric data matrix that can be used in
#' C function. Factors are replaced with their internal values,
#' ordered factors by their internal values or ranks as requested,
#' binary variables are guaranteed to have a zero level. If there are
#' any variables of \code{vtypes = "OP"}, a new column of number of
#' ties for each observation (including tie to itself) is added. Two
#' attributes are added. Attribute \code{"vtypes"} gives the
#' interpreted types as defined above and \code{vtypes = "T"} is used
#' for possible new columns of ties. Attribute \code{scale} gives the
#' scale used in C code. Usually this is the range of variables (and
#' will be ignored for \code{vtypes = "N"}), but for \code{vtypes =
#' "OP"} the scale is adjusted for the number of ties.
#'
#' 
#' @param x Input data frame possibly with different data types
#' @param ordered Policy of interpreting ordered factors (see Details).
#' @param continuous Policy of handling continuous variables (see Details).
#' @param vtypes Variable types for each column of input data. If
#' given, will take precedence over policies in \code{ordered} and
#' \code{continuous}. See Details.
#' 
#' @export
`gowerdistance` <-
    function(x, ordered = "internal", continuous = "metric", vtypes)
{
    if (!is.data.frame(x))
        x <- as.data.frame(x)
    ## allowed variable types ordered as in the C file
    typnam <- c("B","P","N","OP","OR","O","C","T")
    types <- seq_along(typnam)
    names(types) <- typnam
    ## interpret policies
    ordered <- match.arg(ordered,
                         c("internal", "metric", "podani", "nominal"))
    continuous <- match.arg(continuous,
                            c("metric", "presence", "binary"))
    ## Get variable types
    intypes <- sapply(x, function (z) class(z)[1])
    if (!missing(vtypes)) {
        if (length(vtypes) != ncol(x))
            stop("'vtypes' must contain type of every column of 'x'")
    } 
    ## interpret types
    else {
        vtypes <- rep(NA, ncol(x))
        if (any(typ <- intypes == "factor"))
            vtypes[typ] <- "N"
        if (any(typ <- intypes == "ordered"))
            vtypes[typ] <-
                switch(ordered,
                       "internal" = "O",
                       "metric" = "OR",
                       "podani" = "OP",
                       "nominal" = "N")
        ## interpretation of "numeric" can vary with data values
        if (any(typ <- intypes %in% c("numeric", "integer", "logical"))) {
            ## all are treated as metric
            if (continuous == "metric")
                vtypes[typ] <- "C"
            ## non-negative data ignore double-zeros
            if (continuous == "presence") {
                for(i in which(typ))
                    if (min(x[,i] < 0))
                        vtypes[i] <- "C"
                    else
                        vtypes[i] <- "P"
            }
            ## binary data ignore double-zeros
            if (continuous == "binary") {
                for (i in which(typ))
                    if(length(table(x[,i])) == 2)
                        vtypes[i] <- "B"
                    else
                        vtypes[i] <- "C"
            }
        }
    }
    ## start handling data
    x <- data.matrix(x)
    ## replace internal coding or ordered factors with ranks
    if (any(typ <- vtypes %in% c("OR","OP")))
        for(i in which(typ))
            x[,i] <- rank(x[,i])
    ## check that binary variables have minimum=0
    if (any(typ <- vtypes == "B"))
        for(i in which(typ))
            if(min(x[,i] != 0))
                x[,i] <- x[,i] - min(x[,i])
    ## scales for variables are usually their ranges
    scale <- apply(x, 2, function(z) diff(range(z)))
    ## but for Podani scaling we need scale adjustment
    ## -(Timin-1)/2-(Timax-1)/2 that is easiest to find from extreme
    ## averaged ranks
    if (any(typ <- vtypes == "OP"))
        for(i in which(typ))
            scale[i] <- scale[i] - min(x[,i]) + 1 - nrow(x) + max(x[,i])
    ## In Podani scaling we also need the number of tied values for
    ## each observation
    if (any(typ)) {
        ties <- apply(x[, typ, drop = FALSE], 2,
                      function(z) colSums(outer(z, z, "==")))
        x <- cbind(x, ties)
        scale <- c(scale, scale[typ])
        vtypes <- c(vtypes, rep("T", sum(typ)))
    }
    ## return x with attributes we need for Gower distance
    attr(x, "vtypes") <- types[vtypes]
    attr(x, "scale") <- scale
    x
}
