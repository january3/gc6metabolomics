.manu.env <- new.env(parent=emptyenv())
.manu.env$counters <- list()

#' filter out the comments
#' @param lines a character vector
#' @param m_comments pattern for detecting comments; line with this pattern
#'        will be removed
filter_comments <- function(lines, m_comments="^#' -- ") {
  lines[!grepl(m_comments, lines)]
}


#' Insert includes
#' @param lines a character vector
#' @param m_includes pattern to detect include lines
filter_includes <- function(lines, m_includes="^#' @include (.*)", path=NULL) {

  sel <- grep(m_includes, lines)
  fnames <- sapply(lines[sel], function(x) gsub(m_includes, "\\1", x), simplify=FALSE)
  if(!is.null(path)) fnames <- sapply(fnames, function(file) file.path(path, file))

  includes <- sapply(fnames, function(x) readLines(x), simplify=FALSE)

  froms <- c(1, sel + 1) 
  tos   <- c(sel - 1, length(lines))

  lines.2 <- sapply(1:length(froms), function(i) {
    # use unlist(include[i]) rather than includes[[i]], because the last
    # element is outside of the list
    c(lines[froms[i]:tos[i]], unlist(includes[i]))
  }, simplify=FALSE)
  unlist(lines.2)

}


#' Read a text file
#' @param file file name
#' @param path path to file
manu_read <- function(file, path=NULL) {
  if(!is.null(path)) file <- file.path(path, file)
  readLines(file)
}


#' Save lines to a text file
#' @param lines a character vector
#' @param file file name
#' @param path path to the file
#' @return Returns invisibly the path to the file written.
manu_save <- function(lines, file, path=NULL) {
  if(!is.null(path)) file <- file.path(path, file)
  printf("writing to %s", file)
  write(lines, file=file)
  return(invisible(file))
}


#' Process a manuscript
#' @param infile input file
#' @param outfile output file
#' @param path path to files: input, output and include files
#' @param render Whether or not to render (knit) the document immediately
#' @param root.dir Where to process the document. By default, the current directory
#' @return Returns invisibly the path to the manuscript written.
#' @import rmarkdown
manu_process <- function(infile, outfile, 
                         render=FALSE, 
                         output.dir=NULL,
                         knitr.root.dir=NULL, path=NULL
                         ) {
  require(rmarkdown)

  if(is.null(knitr.root.dir)) root.dir <- getwd()
  if(is.null(output.dir)) output.dir <- path

  lines <- manu_read(infile, path)
  # remove comments
  lines <- filter_comments(lines)
  # insert includes
  lines <- filter_includes(lines, path=path)
  outfile <- manu_save(lines, outfile, path)

  if(render) {
    outfile <- render(outfile, knit_root_dir=knitr.root.dir)
  }

  return(invisible(outfile))
}






