library(xml2)
library(rmarkdown)

file = 'inst/2015-08-26-chapter-3-draft.docx'

tmp <- tempfile()
utils::unzip(file, exdir = tmp)  # Unzip to temporary directory
xmlfile <- file.path(tmp, "word", "document.xml")  # Path to xml document
doc_xml = read_xml(xmlfile)
unlink(tmp, recursive = TRUE)  # Delete unzipped files; no longer needed
pars = xml_find_all(doc_xml, "//w:p", xml_ns(doc_xml))
text = xml_text(pars)
tmp2 = tempfile()
cat(text, sep="\n", file = tmp2)
system(paste("fold -w 80 -s", tmp2, " > inst/2015-08-26-chapter-3-draft.rmd"))


