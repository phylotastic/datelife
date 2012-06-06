run <- function(MyTable) {
  # HTML
  out("<!doctype html><html lang='en'><head><title>Data Frame</title>",
      "<meta charset='UTF-8'/></head><body><header><h1><code>textarea</code> ",
      "to <code>data.frame</code></h1></header>", sep = "", eol = "")

  # Retrieve the data in a 'multipart/form-data' form
	oprint(request$c.type)

  out("</body></html>")
  done()
}
