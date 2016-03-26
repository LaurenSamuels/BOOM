PlotQQBOOM <- function(boom, pt.size= 0.5) {
    # boom is an object returned by BOOM()
    # TODO: simplify this after I set other returned lm things to NA if no model
    
    TE.types <- c('est.TEs', 'est.TEs.lm')
    plotlist <- lapply(TE.types, function(tetype){
        if (all(is.na(boom[[tetype]]))) return(NULL)
        dat <- data.frame(y= boom[[tetype]])
        p <- ggplot(dat, aes(sample = y)) +
            geom_point(stat = "qq", size= pt.size,
                dparams= list(mean= mean(dat$y, na.rm= TRUE), sd= sd(dat$y, na.rm= TRUE))) +
            geom_abline(slope= 1, intercept= 0) +
            ylab(tetype)
        p
    })
    if(is.null(plotlist[[2]])) {
        print(plotlist[[1]])
    } else {
        # from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
       # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(1, 2)))
        
        for (i in 1:2) {
            print(plotlist[[i]], vp = viewport(layout.pos.row = 1,
                layout.pos.col = i))
        } 
    }
}
