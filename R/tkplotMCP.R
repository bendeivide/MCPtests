tkplotMCP <- function(results, ...) {
  # Main Window
  base <- tktoplevel(padx=10, pady=10)
  tkwm.geometry(base, "600x400")
  # Change title using tk window manager
  tkwm.title(base,'Plot')
  # Main Window Frame
  main_frame <- tkframe(base, relief="sunken", borderwidth = 1)

  # Canvas
  canvas <- tkcanvas(main_frame)
  # Building
  tkpack(main_frame, canvas, fill = "both", expand=TRUE)


  drawGraph <- function(...) {
    # Graph geometry
    height <- as.numeric(tclvalue(tkwinfo("height", canvas)))
    width <- as.numeric(tclvalue(tkwinfo("width", canvas)))

    # Temp file
    fp <- tempfile(pattern = "MCPtests.",
                   tmpdir = tempdir(),
                   fileext = ".png")

    # Create Image
    png(
      filename = fp,
      width = width,
      height = height,
      units = "px"
    )

    #MCPbarplot(results, ...)
    MCPtests::MCPbarplot(results, ...)

    dev.off()

    tkimage.create("photo", "::image::imgteste", file = fp)

    # Set image to element
    #tkpack.forget(canvas)

    # Create a new canvas with the new image
    #new_canvas <- tkcanvas(main_frame)
    tkcreate(canvas, "image", 0, 0, anchor = "nw", image = "::image::imgteste")

    # Pack new canvas and copy ID to global variable
    # tkpack(new_canvas, fill = "both", expand = TRUE)
    # canvas <- new_canvas
  }


  # CALLBACKS
  onResize <- function() {
    drawGraph()
  }

  # Window Event Binds
  tkbind(base, '<Configure>', onResize)

  tkwait.window(main_frame)
}
