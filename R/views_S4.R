## View manipulation functions

#TODO docs

# --- MistyView class definition
setClass("MistyView", slots = list(name="character", 
                                   abbreviation="character",
                                   building.view.function="closure",
                                   processed="logical"
                                   view.info="list",
                                   view.data="data.frame", 
                                   view.cache="character"))

setGeneric("CallBuildingViewFunction", function(misty.view, force.processing, ...) {
  standardGeneric("CallBuildingViewFunction")
})

setMethod("CallBuildingViewFunction", c(misty.view="MistyView", force.processing...), 
          function(misty.view, force.processing=FALSE, ...) {
          if (!is.null(misty.view@building.view.function) & (!processed | force.processing)) {
              misty.view@building.view.function(...)  
              misty.view@processed=TRUE
            } else {
              message("Processing function for this view is not specified")
            }
            return(misty.view)
          })

# --- MistyPipeline class definition
setClass("MistyPipeline", slots=list(name="character", unique.id="character",
                                      data="data.frame",
                                      cache.location="character", 
                                      views="list"))

setGeneric("InitiatePipeline", function(misty.pipeline, data, verbose=NA, unique.id=NA) {
  standardGeneric("InitiatePipeline")
})

setMethod("InitiatePipeline", c(misty.pipeline = "MistyPipeline", data="data.frame", 
                                verbose="logical", unique.id="character"),
          
          function(misty.pipeline, data, verbose=TRUE, unique.id = NULL) {
            message("Initiating a pipeline")
            
            misty.uniqueid <- ifelse(is.null(unique.id),
                                     GetUniqueId(table),
                                     unique.id)
            
            # create cache
            cache.location <- file.path(".misty.temp", misty.uniqueid)
            CheckAndCreateDirectory(cache.location)
            
            misty.pipeline@unique.id <- misty.uniqueid
            misty.pipeline@data <- data
            misty.pipeline@cache.location <- cache.location
            
            return(misty.pipeline)
          })


setGeneric("AddView", function(misty.pipeline, name, abbreviation,
                               delayed_execution, make_cache, processing_function=NA) {
  standardGeneric("AddView")
})

#TODO learn the changaeble list of params, also data might be something else or even NA? 
setMethod("AddView", c(misty.pipeline="MistyPipeline", 
                       name="character", abbreviation, 
                       delayed_execution="logic",
                       make_cache=TRUE,
                       view.building.funcion="closure", ...),
          #TODO the default function won't work for now
          function(misty.pipeline, data, name, abbreviation,
                   delayed_execution=FALSE, make_cache=TRUE, 
                   view.building.funcion=function(){return(misty.pipeline@data)}, ...) {
            
            message("Adding a new view")
            
            view <- new("MistyView", name=name,  
                        abbreviation=abbreviation)
            
            #Always trying to read from cache first.
            view@data = ReadFromCache(misty.pipeline@cache.location, view@name, view@info)

            if ( !delayed_execution & is.null(view@data)) {
              data <- view.building.funcion(...)
              
              if( is.null(data) )  {
                message("The view is malformed: data is missing. Check the processing function.")
              } else {
                view@processed = TRUE
                view@data = data
                if ( make_cache ) {
                  view@cache.location = WriteToCache(misty.pipeline@cache.location, view@data, view@name, view@info)                
                } 
              }
            
            } else {
              message("Execution of processing on the view is delayed.")
            }
            
            #TODO 
            #msg <- "The new view should have the same number of rows as the intracellular view."
            
            misty.pipeline@views <- c(misty.pipeline@views, name=view)
            return(misty.pipeline)
          })


setGeneric("RemoveViewsByNames", function(misty.pipeline, view.names) {
  standardGeneric("RemoveViewsByNames")
})

setMethod("RemoveViewsByNames", c(misty.pipeline="MistyPipeline", view.names), 
          function(view.names) {
            message("Removing views")
            #TODO this won't work for now
            to.match <- !(view.names %in% c("intracellular", "misty.uniqueid"))
            view.indexes <- match(view.names[to.match], misty.pipeline@views)
            misty.pipeline@views %>% rlist::list.remove(view.indexes)
          })


# ---- View calculating functions 

BuildJuxtaview <- function(expression, positions, verbose=FALSE) {
  # from a deldir object
  get_neighbors <- function(ddobj, id) {
    dplyr::union(
      ddobj$delsgs$ind1[which(ddobj$delsgs$ind2 == id)],
      ddobj$delsgs$ind2[which(ddobj$delsgs$ind1 == id)]
    )
  }
  
  if (verbose) message("Computing triangulation")
  delaunay <- deldir::deldir(as.data.frame(positions))
  
  if (verbose) message("Generating juxtaview")
  
  juxtaview.data <- seq(nrow(expression)) %>% furrr::future_map_dfr(function(cid) {
    alln <- get_neighbors(delaunay, cid)
    # suboptimal placement of dists, but makes conflict if out of scope
    # probably due to lazy evaluations
    dists <- distances::distances(as.data.frame(positions))
    actualn <- alln[which(dists[alln, cid] <= neighbor.thr)]
    data.frame(t(colSums(expression[actualn, ])))
  }, .progress = verbose)
  
  return(juxtaview.data)
}


BuildParaview <- function(expression, verbose=FALSE) { 
  # K.approx is a list containing C, W.plus and s (indexes of sampled columns)
  sample_nystrom_row <- function(K.approx, k) {
    
    # transform k into the row index of reordered K.approx
    k.ind <- which(K.approx$s == k)
    if (purrr::is_empty(k.ind)) {
      k.ind <- length(K.approx$s) + k
    }
    
    cw <- seq(ncol(K.approx$W.plus)) %>%
      purrr::map_dbl(~ K.approx$C[k.ind, ] %*% K.approx$W.plus[, .x])
    
    cwct <- seq(ncol(t(K.approx$C))) %>%
      purrr::map_dbl(~ cw %*% t(K.approx$C)[, .x])
    
    # reorder the columns of cwct so that they correspond to the original order
    cwct[c(K.approx$s, seq_along(cwct)[-s])]
  }
  
  dists <- distances::distances(as.data.frame(positions))
  
  if (is.null(ncells)) {
    if (approx == 1) {
      if (verbose) message("Generating paraview")
      para.view <- seq(nrow(expr)) %>%
        furrr::future_map_dfr(~ data.frame(t(colSums(expr[-.x, ] *
                                                       exp(-(dists[, .x][-.x]^2) / l)))), 
                              .options = furrr::future_options(packages = "distances"))
    }
    else {
      if (approx < 1) approx <- base::round(approx * ncol(dists))
      
      if (verbose) message("Approximating RBF matrix using the Nystrom method")
      # single Nystrom approximation expert, given RBF with paramter l
      s <- sort(sample.int(n = ncol(dists), size = approx))
      C <- exp(-(dists[, s]^2) / l)
      
      # pseudo inverse of W
      W.plus <- MASS::ginv(C[s, ])
      # return Nystrom list
      K.approx <- list(s = s, C = C, W.plus = W.plus)
      
      if (verbose) message("Generating paraview")
      para.view <- seq(nrow(expr)) %>%
        furrr:future_map_dfr(~ data.frame(t(colSums(expr[-.x, ] * sample_nystrom_row(K.approx, .x)[-.x]))))
    }
  } else {
    message("Generating paraview using ", ncells, " nearest neighbors per unit")
    paraview.data <- seq(nrow(expr)) %>%
      furrr::future_map_dfr(function(rowid) {
        knn <- distances::nearest_neighbor_search(dists, ncells + 1, query_indices = rowid)[-1, 1]
        data.frame(t(colSums(expression[knn, ] * exp(-(dists[knn, rowid]^2) / l))),)
      }, .options = furrr::future_options(packages = "distances"))
  }
  
  return(paraview.data)
}

####### Helpers functions

GetUniqueId <- function(object) {
  return(  digest::digest(object, "md5") )
}

CheckAndCreateDirectory <- function(path) {
  if ( !dir.exists(path)  ) {
    dir.create(path, recursive = TRUE, showWarnings = TRUE)
  }
}

ReadFromCache <- function(cache.location, view.name, view.info) {
  filename <- paste0(view.name, view.info, ".rds")
  full.cache.location <- file.path(cache.location, filename) 
  
  if( file.exists(full.cache.location) ) { 
    view.cache <- readr::read_rds(full.cache.location)
  } else {
    message("No cache exists.")
    view.cache <- NULL
  }
  return(view.cache)
}

WriteToCache <- function(cache.location, data_to_write, view.name, view.info="") {
  filename <- paste0(view.name, view.info, ".rds")
  full.cache.location <- file.path(cache.location, filename) 
  
  if ( file.exists(full.cache.location) ) {
    warning("Warning! Rewriting existing cache.")
  } 
  readr::write_rds(data_to_write, full.cache.location)
  return(full.cache.location)
}

######################

RunDefaultPipeline <- function(name="misty_pipeline", data) {
  misty_pipeline <- new("MistyPipeline", name="misty_pipeline")
  misty_pipeline <- InitiatePipeline(misty_pipeline, data)
  
  misty_pipeline <- AddView(misty_pipeline, name="initial_view")
  #TODO specify the params transfered to the functions 
  misty_pipeline <- AddView(misty_pipeline, name="juxtaview", abbreviation="jx", view.building.funcion=BuildJuxtaview)
  misty_pipeline <- AddView(misty_pipeline, name="paraview", abbreviation="para", view.building.funcion=BuildParaview)
  
  return(misty_pipeline)
}
