PlotWeights <- function(boom, xvar= "att.wts.orig", 
    yvar= "BOOM.wts", colorvec= NULL) {
    # boom is an object returned by BOOM()
    dat <- data.frame(x= boom[[xvar]], y= boom[[yvar]],
        # todo: test on things that are already factors
        color.var= if (is.null(colorvec)) NA else factor(colorvec), 
        treat= factor(paste0("Group ", boom$treat)))

    # todo: make sure all same length. or could have
    #   BOOM return the original dset & get the facet var from there

    nicelabels <- c(
        att.wts.orig= "ATT weights from original sample",
        match.wts.orig= "Matching weights from original sample",
        BOOM.wts = "BOOM weights"
    )

    if (is.null(colorvec)) {
        p <- ggplot(data= dat,
            mapping= aes(x=x, y=y)) 
    } else {
        p <- ggplot(data= dat,
            mapping= aes(x=x, y=y, colour= color.var))
        # todo: specify color scale
    }
    p <- p  +
        geom_point(alpha= 0.7) +
        geom_abline(slope= 1, intercept= 0, linetype= "dotted") +
        xlab(nicelabels[xvar]) +
        ylab(nicelabels[yvar]) +
        theme(legend.title=element_blank()) +
        facet_wrap(~treat, ncol=2, scales= "free" )


    print(p)
}
