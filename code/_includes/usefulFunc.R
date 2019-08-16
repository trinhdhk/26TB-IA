#Print a list or a vector as a string

stringify <- function(x, sep, quote = NULL, ...){
  UseMethod('stringify')
}

stringify.default <- function(x, sep = ', ', quote = NULL){
  if (length(quote) > 2) warning('Quotation marks set should <= 2. Only the first 2 is used.')
  if (length(quote) == 2) x <-paste0(quote[1], x, quote[2]) else x <- paste0(quote, x, quote)
  do.call(paste,
          as.list(c(x, sep = sep)))
}

stringify.list <- function(x, sep = c('; ',', '), quote = NULL, setQuote = NULL){
  if (length(sep) == 1) sep <- rep(sep, 2)
  if (length(sep) == 0) sep <- c('; ',', ')
  if (length(setQuote) > 2) warning('Set-quotation-marks set should <= 2. Only the first 2 is used.')
  if (length(setQuote) == 1) setQuote <- rep(setQuote, 2)
  inter.sep <- sep[[1]]
  intra.sep <- sep[[2]]
  atomic.paste <- lapply(x,
                         function(each.x) stringify(each.x, sep = intra.sep, quote = quote))
  atomic.paste <- lapply(atomic.paste,
                         function(each.ap) paste0(setQuote[1], each.ap, setQuote[2]))
  atomic.paste$sep <- inter.sep
  do.call(paste,
          atomic.paste)
}

# Get number of level of a list
n_levels <- function(x){
  n <- 1
  if(any(sapply(x, is.list))){
    n <- n + max(sapply(x, n_levels))
  }
  return(n)
}

# Generate a list of names in list x
names_list <- function(x){
  names <- lapply(seq_along(x), function(i){
    x.1 <- x[i]
    if (is.list(x.1[[1]])) {
      out <- list()
      out[[1]] <-  names(x.1)
      out[[2]] <-  name.list(x.1[[1]])
    } else {
      out <- list()
      out[[1]] <-  names(x.1)
      out[[2]] <- names(x[[i]])
    }
    return(out)
  })
  class(names) <- c('list', 'namesList')
  return(names)
}

# Get the value list of specified level in list x
get_level <- function(x, lv = 1, simplify = TRUE){
  UseMethod('get_level')
}
get_level.list <- function(list, lv = 1, simplify = TRUE){
  out <- sapply(list, function(this){
    
    if (lv == 1) this[[1]]
    else if(is.list(this)) {
      print(this[-1])
      get_level.list(this[-1][[1]], lv-1, simplify = simplify)
    } 
  })
  if (simplify) out <- unlist(out)
  return(out)
}

get_level.namesList <- function(namesList, lv, simplify = TRUE){
  get_level.list(namesList, lv =1, simplify)
}

# Convert a string to an expressiosn
str_to_expr <- function(str){
  pre_parsed <- strsplit(str, split = "\\s?[;,]\\s?", perl = TRUE)
  parsed <- sapply(pre_parsed, function(text) parse(text = text))
  return(parsed)
}

# Generate a list using nesting syntax
quo_to_list <- function(..., eval = TRUE){
  if (!require(rlang)) stop('rlang is needed')
  x.quos <- rlang::enquos(...)
  x <- lapply(x.quos, function(x.quo){
    x.quo <- rlang::quo_squash(x.quo)
    out <- sapply(x.quo[-1], function(q){
      if (is.language(q) & !is.symbol(q)) quo_to_list(UQ(q), eval = eval)
      else {
        o <- if (eval) eval(q) else q
        names(o) <- names(q)
        return(o)
      }
    })
    return(out)
  })
  
  names(x) <- sapply(x.quos, function(x.quo) rlang::quo_squash(x.quo)[[1]])
  return(x)
}
qlist <- function(..., eval = TRUE){
  quo_to_list(..., eval = eval)
}

# Update a list
update_list <- function(x, y, add = TRUE, null.omit = TRUE){
  
  #### Get name list ####
  y.namesList <- names_list(y)
  print(y.namesList)
  
  #### Do the update ####
  for(y.namesList.this in y.namesList){
    if (is.list(y.namesList.this)){
      if (all(is.null(x[y.namesList.this[[1]]][[1]]), add)){
        if (!all(is.null(y[[y.namesList.this[[1]]]])) | !null.omit)
          x[y.namesList.this[[1]]] <-  y[y.namesList.this[[1]]]
      } else { 
        x[[y.namesList.this[[1]]]] <- update_list(x[[y.namesList.this[[1]]]], y[[y.namesList.this[[1]]]], add = add, null.omit = null.omit)
      }
    } else if (add | !is.null(x[[y.namesList.this[[1]]]])) {
      if(!null.omit | !is.null(y[[y.namesList.this]])) x[y.namesList.this] <- y[y.namesList.this] 
    }
  }
  return(x)
}

`update_list<-` <- function(old, add = TRUE, null.omit = TRUE, value)
  update_list(old, value, add = add, null.omit = null.omit)


# Read both csv and xls(x)
read_csv_excel <- function(filePath, ...){
  requireNamespace('readr')
  requireNamespace('readxl')
  ext <- tools::file_ext(filePath)
  if (ext == 'csv') readr::read_csv(filePath, ...)
  else if (ext %in% c('xls', 'xlsx')) readxl::read_excel(filePath, ...)
  else stop('File ext not supported!')
}
