#' Graphical User Interface for MCPtests package
#'
#' \code{guiMCP} A Graphical User Interface (GUI) for
#'     the MCP package
#' @param gui Logical argument, \code{TRUE} or \code{FALSE}. The default is \code{TRUE}
#' @return \code{guiMCP} presents a GUI for the results of the multiple
#'       comparison procedures in the literature. In addition, the GUI
#'       returns a graph of the results, as well as the export of these results
#'       to three types of file extension and latex code.
#' @examples
#' # Loading package
#' library(MCPtests)
#' if (interactive()) {
#'   guiMCP(gui = FALSE)
#' }
#'
#' @import "tcltk"
#' @import "tkrplot"
#' @importFrom "grDevices" "dev.new"
# @importFrom "stats" "aov"
#' @export
guiMCP <- function(gui = TRUE) {
  # Language
  # # portugues
  # Sys.setenv(LANG = "pt_BR")
  #
  # # Ingles
  # Sys.setenv(LANG = "en")

  # Environment of package
  envMCP <- new.env(parent = base::emptyenv())
  assign("dat", NULL, envir = envMCP)
  assign("results", NULL, envir = envMCP)
  assign("LaTeX", NULL, envir = envMCP)
  assign("symmb_promp", NULL, envir = envMCP)

  # Graphical User Interface
  if (gui == TRUE) {

    # Insert Icons
    img <- tclVar(); img2 <- tclVar(); img3 <- tclVar()
    imageinfo <- tkimage.create("photo", img, file = system.file("etc", "info2.gif", package = "MCPtests"))
    imageopen <- tkimage.create("photo", img2, file = system.file("etc", "open.gif", package = "MCPtests"))
    imagesaveas <- tkimage.create("photo", img3, file = system.file("etc", "save-as.gif", package = "MCPtests"))
    ##
    tkimage.create("photo", "::image::logoMCP", file = system.file("etc", "MCP1.gif", package = "MCPtests"))
    tkimage.create("photo", "::image::MCP2", file = system.file("etc", "MCP2.gif", package = "MCPtests"))
    tkimage.create("photo", "::image::MCP3",
                   file = system.file("etc", "MCP3.gif", package = "MCPtests"))
    tkimage.create("photo", "::image::edit",
                   file = system.file("etc", "edit.gif", package = "MCPtests"))
    tkimage.create("photo", "::image::directory",
                   file = system.file("etc", "directory.gif", package = "MCPtests"))


    # Auxiliar functions
    addScrollbars <- function(parent, widget,type=c("both", "x", "y")) {
      if (any(type %in% c("both","x"))) {
        xscr <- ttkscrollbar(parent, orient = "horizontal",
                             command = function(...) tkxview(widget, ...))
        tkconfigure(widget,
                    xscrollcommand = function(...) tkset(xscr,...))
      }

      if(any(type %in% c("both","y"))) {
        yscr <- ttkscrollbar(parent, orient = "vertical",
                             command = function(...) tkyview(widget, ...))
        tkconfigure(widget,
                    yscrollcommand = function(...) tkset(yscr,...))
      }

      ## place in grid
      tkgrid(widget, row=0, column = 0, sticky = "news")
      if (any(type %in% c("both", "x"))) {
        tkgrid(xscr, row=1, column = 0, sticky="ew")
        tkgrid.columnconfigure(parent, 0, weight=1)
      }
      if (any(type %in% c("both", "y"))) {
        tkgrid(yscr,row=0,column=1, sticky="ns")
        tkgrid.rowconfigure(parent, 0, weight=1)
      }


    }

    # ##########################
    # # Configuration of widgets
    # ##########################
    # #Clear the configurations of the option function
    .Tcl("option clear")
    #
    # Button and TButton configurations
    .Tcl("option add *Button.Pady 2
          #option add *Button.Background #e1eff7
          #option add Button.Foreground #e1eff7
          #option add *Button.Foreground black
          option add *Button.Cursor draft_small 46
          option add *TButton.Cursor draft_small 46")
    #
    # # Label configurations
    # .Tcl("option add *Label.Background #e1eff7")
    #
    # # Chackbutton configurations
    # .Tcl("option add *Checkbutton.Background #e1eff7")
    #
    # # Frame configurations
    # .Tcl("option add *Frame.Background #e1eff7")
    #
    #
    # Style TFrame
    #tcl("ttk::style" , "configure" , "Toolbar.TFrame", relief = "solid", expand = TRUE)
    .Tcl("ttk::style configure Toolbar.TFrame -relief solid")
    #
    # # Style LabelFrame
    # .Tcl("ttk::style configure Toolbar.TLabelframe -background #e1eff7
    #       ttk::style configure Toolbar.TLabelframe.Label -background #e1eff7")
    #
    #
    # # Style PanedWindow
    # .Tcl("ttk::style configure Toolbar.TPanedwindow -background #e1eff7")
    #
    # # Class disabled/Enabled
    .Tcl("option add *Ativado.Entry.state normal 81
          option add *Ativado.Label.state normal 81
          option add *Ativado.Button.state normal 81
          option add *Desativado.Entry.state disabled 81
          option add *Desativado.Label.state disabled 81
          option add *Desativado.Button.state disabled 81")


    # Disabled GUI (Type I)
    oldmode <- tclServiceMode(FALSE)

    # Top-level window
    wid <- 1019
    hei <- 700

    topwinMCP <- tktoplevel(
      width = wid,
      height = hei
      #background = "blue"
    )
    # Disabled GUI (Type II)
    #tkwm.state(topwinMCP, "withdraw")

    ####################################
    # Configurations of top-level window
    ####################################

    # Title
    tkwm.title(topwinMCP,
               gettext("GUI to the MCP package", domain = "R-MCPtests"))

    #Icon main toplevel window
    tcl("wm", "iconphoto", topwinMCP, "-default", "::image::logoMCP")

    # Not propagate
    tkpack.propagate(topwinMCP, FALSE)


    # Menu
    # Auxiliar function
    f.read <- NULL
    f.read <- function(file) {
      if (grepl("\\.txt$", file)) {
        if (tclvalue(group_cbox_1_resp) == "TRUE") {
          return(read.table(file, header = TRUE, dec = ",", sep = tclvalue(group_cbox_2_resp)))
        }
        if (tclvalue(group_cbox_1_resp) == "FALSE") {
          return(read.table(file, header = TRUE, sep = tclvalue(group_cbox_2_resp)))
        }
      }
      if (grepl("\\.csv$", file)) {
        if (tclvalue(group_cbox_1_resp)) {
          return(read.table(file, header = TRUE, dec = ",", sep = tclvalue(group_cbox_2_resp)))
        }
        if (tclvalue(group_cbox_1_resp) == FALSE) {
          return(read.table(file, header = TRUE, sep = tclvalue(group_cbox_2_resp)))
        }
      }
    }
    group_cbox_1_resp <- tclVar("FALSE")
    group_cbox_2_resp <- tclVar("")
    #
    openfile <- function(...) {
      tclServiceMode(FALSE)
      if (exists("confdata")) {
        tkwm.deiconify(confdata)
      } else{
        confdata <- tktoplevel()
        tkwm.resizable(confdata, FALSE, FALSE)
        tkwm.title(confdata,
                   gettext("Configurations of the data", domain = "R-MCPtests"))
      }
      # Group of buttons
      tkpack(group_cbox <- tkframe(parent = confdata),
             expand = TRUE, fill = "x", pady = "1m")
      #tkpack.configure(group_cbox, expand = TRUE, fill = "both")

      # Checkbox
      tkpack(group_cbox_1 <- tkcheckbutton(parent = group_cbox,
                                           text = gettext("Comma as decimal points", domain = "R-MCPtests"),
                                           variable = group_cbox_1_resp,
                                           onvalue = "TRUE",
                                           offvalue = "FALSE"),
             anchor = "nw", padx = "1m", side = "left")

      ##
      ##Separator
      tkpack(ttkseparator(parent = group_cbox, orient = "vertical"),
             fill = "both", side = "left")
      ##
      tkpack(tklabel(parent = group_cbox,
                     text = gettext("Separator of variables:",
                                    domain = "R-MCPtests")),
             side = "left", anchor = "nw", padx = "1m"
      )
      ##
      tkpack(group_cbox_2 <- tkentry(textvariable = group_cbox_2_resp,
                                     parent = group_cbox,
                                     width = 5),
             side = "left", anchor = "nw", padx = "1m"
      )

      ##Separator
      tkpack(ttkseparator(parent = confdata, orient = "horizontal"),
             fill = "x")
      tkpack(bconfdata <- ttkbutton(parent = confdata,
                                    text = gettext("Enter the data", domain = "R-MCPtests")), anchor = "e")

      tclServiceMode(TRUE)
      tkfocus(bconfdata)
      funcbconfdata <- function(...){
        filetemp <- tkgetOpenFile(filetypes = paste(
          "{{txt files} {.txt} }" ,
          "{{csv files} {.csv}}" ,
          sep = " "))
        start_dir <- tclvalue(filetemp)
        if (file.exists(start_dir)) {
          tkwm.withdraw(confdata)
          envMCP$dat <- f.read(start_dir)

          tcl(search_results, "delete", "1.0", "end")
          tkinsert(search_results, "end",
                   paste(gettext("The variables in the file:", domain = "R-MCPtests"),
                         "===========================",
                         as.character(start_dir),
                         "===========================",
                         sep = "\n"))
          tkinsert(search_results, "end", "\n")
          tkinsert(search_results, "end",
                   paste(paste(names(f.read(start_dir)), collapse = "\n", sep = "\n"),
                         sep = "\n"), "variablesTag")
          tkinsert(search_results, "end", "\n")
          tkinsert(search_results, "end",
                   paste("==========================="))
          tkmessageBox(message = gettext("Check the data has been loaded correctly. To do this, use the 'Edit/View' button or the 'Output' frame.", domain = "R-MCPtests"))
        }
        if (file.exists(start_dir) == FALSE) {
          tkwm.withdraw(confdata)
          tkmessageBox(message = gettext("No data set has been entered!", domain = "R-MCPtests"))
        }
      }
      tkbind(bconfdata, "<ButtonRelease>", funcbconfdata)
      tkbind(confdata, "<Return>", funcbconfdata)
      tkbind(confdata, "<Escape>", function(){
        tkwm.withdraw(confdata)
      })

    }
    menu_bar <- tkmenu(topwinMCP)
    tkconfigure(topwinMCP, menu = menu_bar)


    # File menu
    file_menu <- tkmenu(menu_bar, tearoff = FALSE)
    chosdir <- function(...) {
      dir_name <- tkchooseDirectory()
      if (nchar(dir_name <- tclvalue(dir_name))) {
        dir_name <- setwd(dir_name)
        on.exit(setwd(dir_name)) # Return initial directory
      }
    }
    tkadd(file_menu, 'command', label = gettext('Choose directory...', domain = "R-MCPtests"),
          accelerator = 'Ctrl+Shift+H', command = chosdir,
          image = "::image::directory", compound = "left")
    tkadd(menu_bar, 'cascade', label = gettext('File', domain = "R-MCPtests"), menu = file_menu)
    tkadd(file_menu, 'command', label = gettext('Open file (.txt or .csv)...', domain = "R-MCPtests"),
          accelerator = 'Ctrl+O', command = openfile, image = tclvalue(imageopen), compound = "left")


    ## Edit menu
    # This variable is important in the event of the "bentry" button
    fedit <- function(...) {
      if (is.null(envMCP$dat)) {
        tkmessageBox(message = gettext("No data set has been entered!", domain = "R-MCPtests"))
      } else{
        envMCP$dat <- edit(envMCP$dat)
      }
    }
    edit_menu <- tkmenu(menu_bar, tearoff = FALSE)
    tkadd(menu_bar, "cascade", label = gettext("Edit", domain = "R-MCPtests"), menu = edit_menu)
    tkadd(edit_menu, "command", label = gettext("Data set...", domain = "R-MCPtests"),
          accelerator = "Ctrl+E", command = fedit,
          image = "::image::edit", compound = "left")


    tkbind(topwinMCP, "<Control-O>", openfile)
    tkbind(topwinMCP, "<Control-o>", openfile)
    tkbind(topwinMCP, "<Control-E>", fedit)
    tkbind(topwinMCP, "<Control-e>", fedit)
    tkbind(topwinMCP, "<Control-Shift-H>", chosdir)
    tkbind(topwinMCP, "<Control-Shift-h>", chosdir)





    ############
    # Overall group
    ############
    #group.all <- NULL
    group.all <- ttkpanedwindow(topwinMCP, orient = "horizontal", style = "Toolbar.TPanedwindow")
    tkpack(group.all, expand = TRUE, fill = "both")
    #tkpack(group.all)
    #tkwinfo("width", group.all)
    #tkwinfo("height", group.all)

    ##############
    # Child groups
    ##############
    #group1 <- NULL
    wid2 <- 400
    group1 <- ttkpanedwindow(group.all, orient = "vertical", width =  wid2, style = "Toolbar.TPanedwindow")
    tkadd(group.all, group1)
    #tkwinfo("width", group1)
    #tkwinfo("height", group1)
    #group2 <- NULL
    group2 <- ttkpanedwindow(group.all, orient = "vertical", width = wid - wid2, style = "Toolbar.TPanedwindow")
    tkadd(group.all, group2)
    #tkwinfo("width", group2)
    #tkwinfo("height", group2)


    # Change width of panes
    #width <- as.integer(tkwinfo("width", group.all))
    #tcl(group.all, "sashpos", 0, floor(0.37 * width))
    #tcl(group.all, "sashpos", 0, 377)


    #######################
    # Child of Child groups
    #######################
    #Q1
    child1.group1 <- ttkframe(parent = group1, style = "Toolbar.TFrame", height = 50,
                              padding = c(3,3,3,3))
    tkadd(group1, child1.group1, weight = 1)

    #Q2
    child2.group1 <- ttkframe(parent = group1, style = "Toolbar.TFrame", height = 200,
                              padding = c(3,3,3,3))
    tkadd(group1, child2.group1, weight = 3)
    #Q3
    # child1.group2 <- ttkframe(parent = group2, style = "Toolbar.TFrame",
    #                           padding = c(3,3,3,3))
    # tkadd(group2, child1.group2)
    #Q4
    child2.group2 <- ttkpanedwindow(parent = group2, orient = "horizontal", width = 669, style = "Toolbar.TPanedwindow")
    tkadd(group2, child2.group2)
    #Q41
    child21.group2 <- ttkframe(parent = child2.group2, style = "Toolbar.TFrame", width = 40,
                               padding = c(3,3,3,3))
    tkadd(child2.group2, child21.group2)
    #Q42
    child22.group2 <- ttkframe(parent = child2.group2, style = "Toolbar.TFrame", width = 300,
                               padding = c(3,3,3,3))
    tkadd(child2.group2, child22.group2)





    ############
    # Working Q1
    ############
    tkpack(Q1 <- tkframe(parent = child1.group1), side = "top", anchor = "n",
           expand = TRUE, fill = "x")

    tkpack(Q1.1 <- tkframe(parent = Q1), side = "top", anchor = "n",
           expand = TRUE, fill = "x")

    #tkpack(tklabel(paren = Q1.1, image = "::image::MCP3"), side = "left", anchor = "nw")
    # Frame
    tkpack(frame.input <- ttklabelframe(text = gettext("Input", domain = "R-MCPtests"), parent = Q1.1,
                                        style = "Toolbar.TLabelframe"),
           side = "left",
           anchor = "n",
           expand = TRUE,
           fill = "x"
    )
    # Group.input.top
    tkpack(group.input.top <- tkframe(parent = frame.input),
           side = "top",
           anchor = "center",
           expand = TRUE,
           fill = "x"
    )
    # Tests
    tkpack(tklabel(text = gettext("Tests:", domain = "R-MCPtests"), parent = group.input.top),
           side = "left",
           anchor = "nw"
    )
    ##
    entry.exp <- c("MGM", "MGR", "SNKM", "TM", "SK")
    vari <- tclVar(entry.exp[1])
    tkpack(tests_box_group <- ttkcombobox(parent = group.input.top,
                                          values = entry.exp,
                                          textvariable = vari,
                                          state = "normal",
                                          justify = "left",
                                          width = 5),
           side = "left",
           anchor = "nw",
           padx = "1m"
    )
    ##
    tkpack(bgit <- tkbutton(text = gettext("help", domain = "R-MCPtests"),
                            parent = group.input.top, borderwidth = 0, image = tclvalue(imageinfo),
                            command = function(...){
                              tkmessageBox(message = gettext("Choose the type of extension for the data output file. If 'latex', the code will be exported to the Console frame. The remaining options will be exported to the selected directory. The choice of the directory can be made in the 'Choosing the directory' button.", domain = "R-MCPtests"))
                            }),
           side = "left",
           anchor = "nw"
    )

    ## Parallelization
    ##
    tkpack(bparalell <- tkbutton(text = gettext("help", domain = "R-MCPtests"),
                                 parent = group.input.top, borderwidth = 0, image = tclvalue(imageinfo),
                                 command = function(...) {
                                   tkmessageBox(message = gettext("Parallelization is still limited for some tests, such as the Scott-Knott's test. However, it is worth mentioning that it is not always efficient to use parallelization, especially when the data set is small.", domain = "R-MCPtests"))
                                 }),
           side = "right",
           anchor = "ne"
    )
    ckeckparallel_resp <- tclVar("FALSE")
    tkpack(ckeckparallel <- tkcheckbutton(parent = group.input.top,
                                          text = gettext("Parallel", domain = "R-MCPtests"),
                                          variable = ckeckparallel_resp,
                                          onvalue = "TRUE",
                                          offvalue = "FALSE"),
           anchor = "ne", padx = "1m", side = "right")

    funpar <- function(...) {
      if (tclvalue(ckeckparallel_resp) == "TRUE") {
        tclvalue(ckeckparallel_resp) <- "FALSE"
      } else {
        tclvalue(ckeckparallel_resp) <- "TRUE"
      }
    }
    tkbind(topwinMCP, "<Alt-P>", funpar)
    tkbind(topwinMCP, "<Alt-p>", funpar)
    ##
    ## Export
    ##
    tkpack(bexp <- tkbutton(text = gettext("help", domain = "R-MCPtests"),
                            parent = group.input.top, borderwidth = 0, image = tclvalue(imageinfo),
                            command = function(...) {
                              tkmessageBox(message = gettext("Choose the type of extension for the data output file. If 'latex', the code will be exported to the Console frame. The remaining options will be exported to the selected directory. The choice of the directory can be made in the 'Choosing the directory' button.", domain = "R-MCPtests"))
                            }),
           side = "right",
           anchor = "ne"
    )
    ##
    ##
    entryexp2 <- c("latex", "txt", "xlsx", "csv")
    vari2 <- tclVar(entryexp2[1])
    tkpack(entry.exp_radio <- ttkcombobox(parent = group.input.top,
                                          values = entryexp2,
                                          textvariable = vari2,
                                          state = "normal",
                                          justify = "left",
                                          width = 8),
           side = "right",
           anchor = "ne",
           padx = "1m"
    )
    ##
    tkpack(tklabel(text = gettext("Extension: ", domain = "R-MCPtests"),
                   parent = group.input.top),
           side = "right",
           anchor = "ne"
    )

    ##Separator
    tkpack(ttkseparator(parent = frame.input, orient = "horizontal"),
           fill = "both", side = "top", anchor = "center", pady = "2m")

    # Data Entry Options
    ####################
    tkpack(group.entry <- tkframe(parent = frame.input),
           side = "top",
           anchor = "center",
           expand = TRUE,
           fill = "x"
    )
    tkpack(tklabel(text = gettext("Options:",
                                  domain = "R-MCPtests"),
                   parent = frame.input),
           side = "left",
           anchor = "nw"
    )
    ##
    entry <- c(gettext("Model", domain = "R-MCPtests"),
               gettext("Response variable", domain = "R-MCPtests"),
               gettext("Averages", domain = "R-MCPtests"))
    vari3 <- tclVar(entry[1])

    # Auxiliar Function

    ##
    tkpack(entry_combobox <- ttkcombobox(parent= frame.input,
                                         values = entry,
                                         textvariable = vari3,
                                         state = "normal",
                                         justify = "left",
                                         width = 15
    ),
    side = "left",
    anchor = "nw",
    padx = "1m"
    )

    tkbind(entry_combobox, "<<ComboboxSelected>>", function(...){
      if (tclvalue(vari3) == gettext("Response variable", domain = "R-MCPtests")) {
        tcl(group.rv, "state", "!disabled")
        child.group.rv <- as.character(tkwinfo("children", group.rv))
        sapply(child.group.rv, function(W) {
          tcl(W, "configure", "-state", "normal")
        })
        child.group.model <- as.character(tkwinfo("children", group.model))
        sapply(child.group.model, function(W) {
          tcl(W, "configure", "-state", "disabled")
        })
        child.group.means <- as.character(tkwinfo("children", groupmeans))
        sapply(child.group.means, function(W) {
          tcl(W, "configure", "-state", "disabled")
        })
      }
      if (tclvalue(vari3) == gettext("Model", domain = "R-MCPtests")) {
        tcl(group.rv, "state", "!disabled")
        child.group.rv <- as.character(tkwinfo("children", group.rv))
        sapply(child.group.rv, function(W) {
          tcl(W, "configure", "-state", "disabled")
        })
        child.group.model <- as.character(tkwinfo("children", group.model))
        sapply(child.group.model, function(W) {
          tcl(W, "configure", "-state", "normal")
        })
        child.group.means <- as.character(tkwinfo("children", groupmeans))
        sapply(child.group.means, function(W) {
          tcl(W, "configure", "-state", "disabled")
        })
      }
      if (tclvalue(vari3) == gettext("Averages", domain = "R-MCPtests")) {
        tcl(group.rv, "state", "!disabled")
        child.group.rv <- as.character(tkwinfo("children", group.rv))
        sapply(child.group.rv, function(W) {
          tcl(W, "configure", "-state", "disabled")
        })
        child.group.model <- as.character(tkwinfo("children", group.model))
        sapply(child.group.model, function(W) {
          tcl(W, "configure", "-state", "disabled")
        })
        child.group.means <- as.character(tkwinfo("children", groupmeans))
        sapply(child.group.means, function(W) {
          tcl(W, "configure", "-state", "normal")
        })
      }
    })

    tkpack(bentry <- tkbutton(text = gettext("help",
                                             domain = "R-MCPtests"),
                              parent = frame.input, borderwidth = 0, image = tclvalue(imageinfo),
                              command = function(...){
                                tkmessageBox(message = gettext("Choose the type of data entry. If 'Model', enter the experimental model. If 'Response Variable', enter the object name of the responses and treatments variables and if 'Averages', enter the vector of means and treatments. For more details, use the help button on each option. \n The 'Averages' option is the only one that will not need to load a data set. Just enter the recommended options and then click the 'Calculate' button.", domain = "R-MCPtests"))
                              }
    ),
    side = "left",
    anchor = "nw"
    )
    ##


    # Significance level
    tkpack(bsl <- tkbutton(text = gettext("help", domain = "R-MCPtests"),
                           parent = frame.input, borderwidth = 0, image = tclvalue(imageinfo),
                           command = function(...){
                             tkmessageBox(message = gettext("Enter the value of the significance level. This value is a number between 0 and 1.", domain = "R-MCPtests"))
                           }),
           side = "right",
           anchor = "ne"
    )
    ##
    ##
    vari4 <- tclVar("0.05")
    ##
    tkpack(slentry <- tkentry(parent = frame.input,
                              textvariable = vari4, width = 5),
           anchor = "ne", padx = "2m", side = "right")

    tkpack(tklabel(text = gettext("Significance:",
                                  domain = "R-MCPtests"),
                   parent = frame.input),
           side = "right",
           anchor = "ne"
    )

    # Output
    tkpack(frame.output <- ttklabelframe(text = gettext("Output:",
                                                        domain = "R-MCPtests"),
                                         parent = Q1,
                                         width = 30, height = 10, style = "Toolbar.TLabelframe"),
           side = "top",
           anchor = "n",
           expand = TRUE,
           fill = "both"
    )

    ##
    tkpack(search_results <- tktext(parent = frame.output, padx = "1m", pady = "1m",
                                    height = 20, wrap = "none"),
           side = "top",
           anchor = "n",
           expand = TRUE,
           fill = "both")
    tktag.configure(search_results, "variablesTag", foreground = "blue" ,
                  font = "courier 10 italic" )

    addScrollbars(frame.output, search_results)


    ############
    # Working Q2
    ############

    # Button Calculate
    ##################
    tkpack(calculate_button <- ttkbutton(text = gettext("Calculate",
                                                        domain = "R-MCPtests"),
                                         parent = child2.group1),
           side = "bottom",
           anchor = "s",
           expand = TRUE,
           fill = "x"
    )
    ##
    # Auxiliar function (Plot)
    # Building function based in the tkrplot package

    tkplot <- function(parent, fun, hscale=1, vscale=1, ...) {
      # tkrplot::.make.tkindex()
      image <- paste("Rplot", .make.tkindex(), sep = "")
      # tkrplot::.my.tkdev()
      .my.tkdev(hscale, vscale)
      try(fun())
      .Tcl(paste("image create Rplot", image))
      lab <- tktext(parent, height = 50, wrap = "none")
      tcl(lab, "image", "create", "0.0", "-image", image, "-align", "center")
      tkbind(lab,"<Destroy>", function() .Tcl(paste("image delete", image)))
      lab$image <- image
      lab$fun <- fun
      lab$hscale <- hscale
      lab$vscale <- vscale
      lab
    }
    ##
    calcular <- function(...){
      tcl(console, "delete", "1.0", "end")
      if (tclvalue(vari3) == gettext("Model", domain = "R-MCPtests")) {
        capture.output(envMCP$results <- MCPtests::MCPtest(y = aov(eval(parse(text = tclvalue(vari5))), data = envMCP$dat),
                                                                  trt = tclvalue(vari6),
                                                                  alpha = eval(parse(text = tclvalue(vari4))),
                                                                  MCP = tclvalue(vari),
                                                                  parallel = eval(parse(text = tclvalue(ckeckparallel_resp)))))
        objtreat <- as.factor(envMCP$dat[,tclvalue(vari6)])

        # Results in the Console
        if (tclvalue(vari2) == "latex") {
          envMCP$LaTeX <- MCPtests::MCPwrite(x = envMCP$results, extension = 'latex')
          eval_cmd_chunk(console, "results; LaTeX")
        }
        if (tclvalue(vari2) != "latex") {
          MCPtests::MCPwrite(x = envMCP$results, MCP = 'all', extension = tclvalue(vari2), dataMR = 'all')
          eval_cmd_chunk(console, "results")
        }

        ##
        #Update plots
        sapply(as.character(tkwinfo("children", frame.graf)),
               function(W) tcl("destroy", W))

        #Frame label
        tkpack(framelabel <- ttkframe(frame.graf, style = "Toolbar.TFrame",
                                      padding = c(3,3,3,3)),
               side = "top", fill = "x", expand = TRUE)

        # Frame label aux
        tkpack(framelabelaux <- ttkframe(framelabel, style = "Toolbar.TFrame",
                                         padding = c(3,3,3,3)),
               side = "top", fill = "x", expand = TRUE)

        #Frame plot
        tkpack(frameplot <- ttkframe(frame.graf, style = "Toolbar.TFrame",
                                     padding = c(3,3,3,3), height = 360),
               side = "top", fill = "both", expand = TRUE)

        # Creating an image in the tcltk
        fator <- 1
        #Plot
        ##
        plottopwinMCP <- tkplot(parent = frameplot,
                                function(...) {
                                  MCPtests::MCPbarplot(envMCP$results,col = heat.colors(length(levels(objtreat))))
                                }, hscale = fator, vscale = fator)
        # Auxiliar function
        f <- function(...) {
          fatorup <- as.numeric(tclvalue("fator"))
          if (fator != fatorup) {
            fator <<- fatorup
            tkrreplot(plottopwinMCP, hscale = fator, vscale = fator)
          }
        }
        # Button save as...
        #graphics.off() # Erasing All Graphics Devices
        dispplot <- NULL
        bsaveas <- tkbutton(parent = framelabelaux, text = gettext("Save as...", domain = "R-MCPtests"),
                            borderwidth = 0, underline = 0,
                            image = tclvalue(imagesaveas), compound = "top",
                            command = function(...){
                              dispplot <<- grDevices::dev.new(noRStudioGD = TRUE) # New device plot
                              #grDevices::dev.new(noRStudioGD = TRUE) # New device plot
                              if (tclvalue(vari3) == gettext("Model", domain = "R-MCPtests")) {
                                objtreat <- as.factor(envMCP$dat[,tclvalue(vari6)])
                              }
                              if (tclvalue(vari3) == gettext("Response variable", domain = "R-MCPtests")) {
                                objtreat <- as.factor(envMCP$dat[,tclvalue(vari8)])
                              }
                              if (tclvalue(vari3) == gettext("Averages", domain = "R-MCPtests")) {
                                # Treatment levels
                                trat <- strsplit(tclvalue(vari12), split = ",", perl = TRUE)[[1]]
                                trat <- as.factor(trat)
                                objtreat <- trat
                              }
                              xlab <- NULL; ylab <- NULL
                              eval(parse(text = tclvalue(vari18)))
                              hor <- if (tclvalue(vari17) == "") FALSE else eval(parse(text = tclvalue(vari17)))
                              if (tclvalue(vari16) == "") {
                                color <- heat.colors(length(levels(objtreat)))
                              } else {
                                if (any(c(grep("\'", tclvalue(vari16)), grep("\"", tclvalue(vari16))) == 1)) {
                                  if (length(grep("\'", tclvalue(vari16))) > 0) {
                                    color <- gsub("\'", "", tclvalue(vari16))
                                  }
                                  if (length(grep("\"", tclvalue(vari16))) > 0) {
                                    color <- gsub("\"", "", tclvalue(vari16))
                                  }
                                } else{
                                  color <- eval(parse(text = tclvalue(vari16)))
                                }
                              }
                              MCPtests::MCPbarplot(envMCP$results, col = color, horiz = hor, xlab = xlab, ylab = ylab)
                            })
        tkpack(bsaveas, side = "left")

        # Frame label aux2
        tkpack(framelabelaux2 <- ttkframe(framelabelaux, style = "Toolbar.TFrame",
                                          padding = c(3,3,3,3)),
               side = "top", fill = "x", expand = TRUE)
        #Label
        scalelabel <- tklabel(parent = framelabelaux2, text = gettext("Scale of Plot", domain = "R-MCPtests"))
        # Scale of plot
        s <- tkscale(framelabelaux2, command = f, from = 1, to = 3.00, variable = "fator",
                     showvalue = TRUE, resolution = 0.05, orient = "horiz", borderwidth = 0)

        # Drawing the widget
        tkpack(s)
        tkpack(scalelabel, before = s)
        tkpack.configure(s, fill = "both", expand = TRUE, side = "bottom", ipady = 5)
        #tkpack.configure(scalelabel, side = "bottom", anchor = "s")



        #tkgrid(plottopwinMCP,  row = 0, column = 0, sticky = "news")
        tkpack(plottopwinMCP,  fill = "both", expand = TRUE)
        addScrollbars(frameplot, plottopwinMCP)

        # #addScrollbars(frame.graf, img)

      }
      if (tclvalue(vari3) == gettext("Response variable", domain = "R-MCPtests")) {
        objrv <- dat[,tclvalue(vari7)]
        objtreat <- as.factor(dat[,tclvalue(vari8)])
        capture.output(results <- results <<- MCPtests::MCPtest(y = objrv,
                                                                trt = objtreat,
                                                                dferror = eval(parse(text = tclvalue(vari9))),
                                                                mserror = eval(parse(text = tclvalue(vari10))),
                                                                alpha = eval(parse(text = tclvalue(vari4))),
                                                                MCP = tclvalue(vari)))

        # Results in the Console
        if (tclvalue(vari2) == "latex") {
          eval_cmd_chunk(console, "results; MCPtests::MCPwrite(x = results, extension = 'latex')")
        }
        if (tclvalue(vari2) != "latex") {
          MCPtests::MCPwrite(x = results, MCP = 'all', extension = tclvalue(vari2), dataMR = 'all')
          eval_cmd_chunk(console, "results")
        }

        ##
        #Update plots
        sapply(as.character(tkwinfo("children", frame.graf)),
               function(W) tcl("destroy", W))

        #Frame label
        tkpack(framelabel <- ttkframe(frame.graf, style = "Toolbar.TFrame",
                                      padding = c(3,3,3,3)),
               side = "top", fill = "x", expand = TRUE)

        # Frame label aux
        tkpack(framelabelaux <- ttkframe(framelabel, style = "Toolbar.TFrame",
                                         padding = c(3,3,3,3)),
               side = "top", fill = "x", expand = TRUE)

        #Frame plot
        tkpack(frameplot <- ttkframe(frame.graf, style = "Toolbar.TFrame",
                                     padding = c(3,3,3,3), height = 360),
               side = "top", fill = "both", expand = TRUE)

        # Creating an image in the tcltk
        fator <- 1
        #Plot
        ##
        plottopwinMCP <- tkplot(parent = frameplot,
                                function(...) {
                                  MCPtests::MCPbarplot(envMCP$results,col = heat.colors(length(levels(objtreat))))
                                }, hscale = fator, vscale = fator)
        # Auxiliar function
        f <- function(...) {
          fatorup <- as.numeric(tclvalue("fator"))
          if (fator != fatorup) {
            fator <<- fatorup
            tkrreplot(plottopwinMCP, hscale = fator, vscale = fator)
          }
        }
        # Button save as...
        #graphics.off() # Erasing All Graphics Devices
        dispplot <- NULL
        bsaveas <- tkbutton(parent = framelabelaux, text = gettext("Save as...", domain = "R-MCPtests"),
                            borderwidth = 0, underline = 0,
                            image = tclvalue(imagesaveas), compound = "top",
                            command = function(...){
                              dispplot <<- grDevices::dev.new(noRStudioGD = TRUE) # New device plot
                              #grDevices::dev.new(noRStudioGD = TRUE) # New device plot
                              if (tclvalue(vari3) == gettext("Model", domain = "R-MCPtests")) {
                                objtreat <- as.factor(envMCP$dat[,tclvalue(vari6)])
                              }
                              if (tclvalue(vari3) == gettext("Response variable", domain = "R-MCPtests")) {
                                objtreat <- as.factor(envMCP$dat[,tclvalue(vari8)])
                              }
                              if (tclvalue(vari3) == gettext("Averages", domain = "R-MCPtests")) {
                                # Treatment levels
                                trat <- strsplit(tclvalue(vari12), split = ",", perl = TRUE)[[1]]
                                trat <- as.factor(trat)
                                objtreat <- trat
                              }
                              xlab <- NULL; ylab <- NULL
                              eval(parse(text = tclvalue(vari18)))
                              hor <- if (tclvalue(vari17) == "") FALSE else eval(parse(text = tclvalue(vari17)))
                              if (tclvalue(vari16) == "") {
                                color <- heat.colors(length(levels(objtreat)))
                              } else {
                                if (any(c(grep("\'", tclvalue(vari16)), grep("\"", tclvalue(vari16))) == 1)) {
                                  if (length(grep("\'", tclvalue(vari16))) > 0) {
                                    color <- gsub("\'", "", tclvalue(vari16))
                                  }
                                  if (length(grep("\"", tclvalue(vari16))) > 0) {
                                    color <- gsub("\"", "", tclvalue(vari16))
                                  }
                                } else{
                                  color <- eval(parse(text = tclvalue(vari16)))
                                }
                              }
                              MCPtests::MCPbarplot(envMCP$results, col = color, horiz = hor, xlab = xlab, ylab = ylab)
                            })
        tkpack(bsaveas, side = "left")

        # Frame label aux2
        tkpack(framelabelaux2 <- ttkframe(framelabelaux, style = "Toolbar.TFrame",
                                          padding = c(3,3,3,3)),
               side = "top", fill = "x", expand = TRUE)
        #Label
        scalelabel <- tklabel(parent = framelabelaux2, text = gettext("Scale of Plot", domain = "R-MCPtests"))
        # Scale of plot
        s <- tkscale(framelabelaux2, command = f, from = 1, to = 3.00, variable = "fator",
                     showvalue = TRUE, resolution = 0.05, orient = "horiz", borderwidth = 0)

        # Drawing the widget
        tkpack(s)
        tkpack(scalelabel, before = s)
        tkpack.configure(s, fill = "both", expand = TRUE, side = "bottom", ipady = 5)
        #tkpack.configure(scalelabel, side = "bottom", anchor = "s")



        #tkgrid(plottopwinMCP,  row = 0, column = 0, sticky = "news")
        tkpack(plottopwinMCP,  fill = "both", expand = TRUE)
        addScrollbars(frameplot, plottopwinMCP)
      }
      if (tclvalue(vari3) == gettext("Averages", domain = "R-MCPtests")) {
        # if (!is.factor(eval(parse(text = svalue(gme2d))))){
        #   gmessage("The trt argument must be factor")
        # }
        # Mean vector
        aver <- paste("c(", tclvalue(vari11), ")")
        # Treatment levels
        trat <- strsplit(tclvalue(vari12), split = ",", perl = TRUE)[[1]]
        trat <- as.factor(trat)

        # Functions
        capture.output(results <- results <<- MCPtests::MCPtest(y = eval(parse(text = aver)),
                                                                trt = trat,
                                                                dferror = eval(parse(text = tclvalue(vari13))),
                                                                mserror = eval(parse(text = tclvalue(vari14))),
                                                                alpha = eval(parse(text = tclvalue(vari4))),
                                                                MCP = tclvalue(vari),
                                                                replication = eval(parse(text = tclvalue(vari15))),
                                                                ismean = TRUE))

        # Results in the Console
        if (tclvalue(vari2) == "latex") {
          eval_cmd_chunk(console, "results; MCPtests::MCPwrite(x = results, extension = 'latex')")
        }
        if (tclvalue(vari2) != "latex") {
          MCPtests::MCPwrite(x = results, MCP = 'all', extension = tclvalue(vari2), dataMR = 'all')
          eval_cmd_chunk(console, "results")
        }

        ##
        #Update plots
        sapply(as.character(tkwinfo("children", frame.graf)),
               function(W) tcl("destroy", W))

        #Frame label
        tkpack(framelabel <- ttkframe(frame.graf, style = "Toolbar.TFrame",
                                      padding = c(3,3,3,3)),
               side = "top", fill = "x", expand = TRUE)

        # Frame label aux
        tkpack(framelabelaux <- ttkframe(framelabel, style = "Toolbar.TFrame",
                                         padding = c(3,3,3,3)),
               side = "top", fill = "x", expand = TRUE)

        #Frame plot
        tkpack(frameplot <- ttkframe(frame.graf, style = "Toolbar.TFrame",
                                     padding = c(3,3,3,3), height = 360),
               side = "top", fill = "both", expand = TRUE)

        # Creating an image in the tcltk
        fator <- 1
        #Plot
        ##
        plottopwinMCP <- tkplot(parent = frameplot,
                                function(...) {
                                  MCPtests::MCPbarplot(envMCP$results,col = heat.colors(length(levels(trat))))
                                }, hscale = fator, vscale = fator)
        # Auxiliar function
        f <- function(...) {
          fatorup <- as.numeric(tclvalue("fator"))
          if (fator != fatorup) {
            fator <<- fatorup
            tkrreplot(plottopwinMCP, hscale = fator, vscale = fator)
          }
        }
        # Button save as...
        #graphics.off() # Erasing All Graphics Devices
        dispplot <- NULL
        bsaveas <- tkbutton(parent = framelabelaux, text = gettext("Save as...", domain = "R-MCPtests"),
                            borderwidth = 0, underline = 0,
                            image = tclvalue(imagesaveas), compound = "top",
                            command = function(...){
                              dispplot <<- grDevices::dev.new(noRStudioGD = TRUE) # New device plot
                              #grDevices::dev.new(noRStudioGD = TRUE) # New device plot
                              if (tclvalue(vari3) == gettext("Model", domain = "R-MCPtests")) {
                                objtreat <- as.factor(envMCP$dat[,tclvalue(vari6)])
                              }
                              if (tclvalue(vari3) == gettext("Response variable", domain = "R-MCPtests")) {
                                objtreat <- as.factor(envMCP$dat[,tclvalue(vari8)])
                              }
                              if (tclvalue(vari3) == gettext("Averages", domain = "R-MCPtests")) {
                                # Treatment levels
                                trat <- strsplit(tclvalue(vari12), split = ",", perl = TRUE)[[1]]
                                trat <- as.factor(trat)
                                objtreat <- trat
                              }
                              xlab <- NULL; ylab <- NULL
                              eval(parse(text = tclvalue(vari18)))
                              hor <- if (tclvalue(vari17) == "") FALSE else eval(parse(text = tclvalue(vari17)))
                              if (tclvalue(vari16) == "") {
                                color <- heat.colors(length(levels(objtreat)))
                              } else {
                                if (any(c(grep("\'", tclvalue(vari16)), grep("\"", tclvalue(vari16))) == 1)) {
                                  if (length(grep("\'", tclvalue(vari16))) > 0) {
                                    color <- gsub("\'", "", tclvalue(vari16))
                                  }
                                  if (length(grep("\"", tclvalue(vari16))) > 0) {
                                    color <- gsub("\"", "", tclvalue(vari16))
                                  }
                                } else{
                                  color <- eval(parse(text = tclvalue(vari16)))
                                }
                              }
                              MCPtests::MCPbarplot(envMCP$results, col = color, horiz = hor, xlab = xlab, ylab = ylab)
                            })
        tkpack(bsaveas, side = "left")

        # Frame label aux2
        tkpack(framelabelaux2 <- ttkframe(framelabelaux, style = "Toolbar.TFrame",
                                          padding = c(3,3,3,3)),
               side = "top", fill = "x", expand = TRUE)
        #Label
        scalelabel <- tklabel(parent = framelabelaux2, text = gettext("Scale of Plot", domain = "R-MCPtests"))
        # Scale of plot
        s <- tkscale(framelabelaux2, command = f, from = 1, to = 3.00, variable = "fator",
                     showvalue = TRUE, resolution = 0.05, orient = "horiz", borderwidth = 0)

        # Drawing the widget
        tkpack(s)
        tkpack(scalelabel, before = s)
        tkpack.configure(s, fill = "both", expand = TRUE, side = "bottom", ipady = 5)
        #tkpack.configure(scalelabel, side = "bottom", anchor = "s")



        #tkgrid(plottopwinMCP,  row = 0, column = 0, sticky = "news")
        tkpack(plottopwinMCP,  fill = "both", expand = TRUE)
        addScrollbars(frameplot, plottopwinMCP)
      }
    }
    tkbind(calculate_button, "<ButtonRelease>", calcular)
    tkbind(topwinMCP, "<Control-Return>", calcular)

    # Console
    #########
    tkpack(frame.console <- ttklabelframe(text = gettext("Console:",
                                                         domain = "R-MCPtests"),
                                          parent = child2.group1, style = "Toolbar.TLabelframe"),
           side = "bottom", expand = TRUE, fill = "x", anchor = "s")
    ##
    ##
    tkpack(console <- tktext(parent = frame.console, height = 40, wrap = "none"),
           side = "bottom", anchor = "s", expand = TRUE, fill = "both")
    ##
    ##
    addScrollbars(frame.console, console)
    ##
    tktag.configure(console, "commandTag", foreground = "blue" ,
                    font = "courier 10 italic" )
    tktag.configure(console, "outputTag", font = "courier 10" )
    tktag.configure(console, "errorTag", foreground = "red" ,
                     font = "courier 10 bold" )
    ##
    # Console function
    eval_cmd_chunk <- function(console , cmds) {
      cmd_chunks <- try(parse(text = cmds), silent = TRUE)
      if (inherits(cmd_chunks , "try-error")) {
        tkinsert(t,"end", "Error", "errorTag") # add markup tag
      }
      for (cmd in cmd_chunks) {
        cutoff <- 0.75 * getOption("width")
        dcmd <- deparse(cmd , width.cutoff = cutoff)
        command <-
          paste("MCPtests> ",
                paste(dcmd, collapse = paste("\n", getOption("continue"), sep = "")) ,
                sep = "" , collapse = "")
        tkinsert(console , "end" , command , "commandTag" )
        tkinsert(console , "end" , "\n")
        ## output, should check for errors in eval!
        # The function sink() stores the output in a file
        output <- capture.output(eval(cmd, envir = envMCP))
        output <- paste(output , collapse = "\n" )
        tkinsert(console, "end", output, "outputTag" )
        tkinsert(console, "end" , "\n")
      }
    }





    ############
    # Working Q3
    ############

    # Logo MCP
    tkpack(tklabel(parent = child21.group2, image = "::image::MCP3"))

    # Option 'Model'
    ################
    tkpack(group.model <- ttklabelframe(text = gettext("Model",
                                                       domain = "R-MCPtests"),
                                        parent = child21.group2, style = "Toolbar.TLabelframe"),
           expand = TRUE, fill = "both")
    ## Enter model
    tkgrid(label.model <- tklabel(text = gettext("Enter model:",
                                  domain = "R-MCPtests"),
                   parent = group.model),
           row = 0, column = 0, sticky = "e"
    )
    ##
    vari5 <- tclVar("")
    tkgrid(gm1d <- tkentry(textvariable = vari5,
                            parent = group.model, width = 20),
           row = 0, column = 1, sticky = "e", padx = 2
    )
    ##
    tkgrid(bm1d <- tkbutton(text = gettext("help", domain = "R-MCPtests"),
                             parent = group.model, borderwidth = 0, image = tclvalue(imageinfo),
                            command = function(...) {
                              tkmessageBox(message = gettext("Enter the experimental model of type Response Variable (RV) ~ Predictive Variables (PV). These variables are in the 'Output' frame, after entering the data set. For example, in a randomized block design, assuming the 'treat' object corresponding to the treatments, 'block' object corresponding to the blocks and 'resp' object corresponding to the variable response. So, you must enter the following expression: resp ~ trat + block.", domain = "R-MCPtests"))
                            }),
           row = 0, column = 2, sticky = "e")
    ##
    ## Treatment
    tkgrid(lm.treat <- tklabel(text = gettext("Treatment:", domain = "R-MCPtests"),
                   parent = group.model),
           row = 1, column = 0, sticky = "e"
    )
    ##
    vari6 <- tclVar("")
    tkgrid(treat <- tkentry(textvariable = vari6,
                             parent = group.model, width = 20),
           row = 1, column = 1, sticky = "e", padx = 2
    )
    ##
    tkgrid(bm2d <- tkbutton(text = gettext("help", domain = "R-MCPtests"),
                             parent = group.model, borderwidth = 0, image = tclvalue(imageinfo),
                            command = function(...) {
                              tkmessageBox(message = gettext("Enter the name of the treatments in the experiment model inserted above, in Predictive Variables (PV). The name of the treatment are in the 'Output' frame, after entering the data set. Inserted all the arguments above, click on the 'Calculate' button.", domain = "R-MCPtests"))
                            }),
           row = 1, column = 2, sticky = "e"
    )
    ##

    # Option 'Response variable'
    ############################
    tkpack(group.rv <- ttklabelframe(text = gettext("Response variable",
                                                    domain = "R-MCPtests"),
                                     parent = child21.group2, style = "Toolbar.TLabelframe",
                                     class = "Desativado"),
           expand = TRUE, fill = "both")
    # Disabled group.rv
    tcl(group.rv, "state", "disabled")
    ## Response
    tkgrid(label.rv <- tklabel(text = gettext("Response:",
                                  domain = "R-MCPtests"),
                   parent = group.rv),
           row = 0, column = 0, sticky = "e"
    )
    ##
    vari7 <- tclVar("")
    tkgrid(grv1d <- tkentry(textvariable = vari7,
                             parent = group.rv, width = 20),
           row = 0, column = 1, sticky = "e", padx = 2
    )
    ##

    tkgrid(bgrv1d <- tkbutton(text = gettext("help", domain = "R-MCPtests"),
                               parent = group.rv, borderwidth = 0, image = tclvalue(imageinfo),
                              command = function(...) tkmessageBox(message = gettext("Insert the name of the variable response of the experiment model. The name of the response variable are in the 'Output' frame, after inserting the data set.", domain = "R-MCPtests"))),
           row = 0, column = 2, sticky = "e")
    ##
    # tkbind(bgrv1d, "<ButtonRelease>", function(){
    #   tkmessageBox(message = gettext("Insert the name of the variable response of the experiment model. The name of the response variable are in the 'Output' frame, after inserting the data set.", domain = "R-MCPtests"))
    # })
    ## Treatment
    tkgrid(lrv.treat <- tklabel(text = gettext("Treatment:", domain = "R-MCPtests"),
                   parent = group.rv),
           row = 1, column = 0, sticky = "e"
    )
    ##
    vari8 <- tclVar("")
    tkgrid(grv2d <- tkentry(textvariable = vari8,
                             parent = group.rv, width = 20),
           row = 1, column = 1, sticky = "e", padx = 2
    )
    ##
    tkgrid(bgrv2d <- tkbutton(text = gettext("help", domain = "R-MCPtests"),
                               parent = group.rv, borderwidth = 0, image = tclvalue(imageinfo),
                              command = function(...) tkmessageBox(message = gettext("Enter the name of the treatment of the experiment model. The name of the response variable are in the 'Output' frame, after entering the data set.", domain = "R-MCPtests"))),
           row = 1, column = 2, sticky = "e"
    )
    ##

    ## DFerror
    tkgrid(lrv.dfe <- tklabel(text = gettext("DFerror:", domain = "R-MCPtests"),
                   parent = group.rv),
           row = 2, column = 0, sticky = "e"
    )
    ##
    vari9 <- tclVar("")
    tkgrid(grv3d <- tkentry(textvariable = vari9,
                             parent = group.rv, width = 20),
           row = 2, column = 1, sticky = "e", padx = 2
    )
    ##
    tkgrid(bgrv3d <- tkbutton(text = gettext("help", domain = "R-MCPtests"),
                               parent = group.rv, borderwidth = 0, image = tclvalue(imageinfo),
                              command = function(...) {
                                tkmessageBox(message = gettext("Enter the numeric value of the degrees of freedom of the mean square error of the experiment model.", domain = "R-MCPtests"))
                              }),
           row = 2, column = 2, sticky = "e"
    )
    ##
    ## MSerror
    tkgrid(lrv.mse <- tklabel(text = gettext("MSerror:", domain = "R-MCPtests"),
                   parent = group.rv),
           row = 3, column = 0, sticky = "e"
    )
    ##
    vari10 <- tclVar("")
    tkgrid(grv4d <- tkentry(textvariable = vari10,
                             parent = group.rv, width = 20),
           row = 3, column = 1, sticky = "e", padx = 2
    )
    ##
    tkgrid(bgrv4d <- tkbutton(text = gettext("help", domain = "R-MCPtests"),
                               parent = group.rv, borderwidth = 0, image = tclvalue(imageinfo),
                              command = function(...) {
                                tkmessageBox(message = gettext("Enter the value of the mean square error of the experiment model. The value is numeric. Inserted all the arguments above, click on the 'Calculate' button.", domain = "R-MCPtests"))
                              }),
           row = 3, column = 2, sticky = "e"
    )
    ##

    # Option 'Averages'
    ############################
    tkpack(groupmeans <- ttklabelframe(text = gettext("Averages",
                                                      domain = "R-MCPtests"),
                                       parent = child21.group2, style = "Toolbar.TLabelframe",
                                       class = "Desativado"),
           expand = TRUE, fill = "both")
    # Disabled groupmeans
    tcl(groupmeans, "state", "disabled")
    ## Averages
    tkgrid(label.averages <- tklabel(text = gettext("Averages:",
                                  domain = "R-MCPtests"),
                   parent = groupmeans),
           row = 0, column = 0, sticky = "e"
    )
    ##
    vari11 <- tclVar("")
    tkgrid(gme1d <- tkentry(textvariable = vari11,
                             parent = groupmeans, width = 20),
           row = 0, column = 1, sticky = "e", padx = 2
    )
    ##
    ##
    tkgrid(bgme1d <- tkbutton(text = gettext("help", domain = "R-MCPtests"),
                               parent = groupmeans, borderwidth = 0, image = tclvalue(imageinfo),
                              command = function(...) {
                                tkmessageBox(message = gettext("Enter the values of the averages. Each mean of the vector must be separated by a comma. For example, for the vector of the average of four treatments: 10, 20, 30, 40. You do not need to use the concatenate function, i.e., c().", domain = "R-MCPtests"))
                              }),
           row = 0, column = 2, sticky = "e"
    )
    ##
    ## Treatment
    tkgrid(la.treat <- tklabel(text = gettext("Treatment:", domain = "R-MCPtests"),
                   parent = groupmeans),
           row = 1, column = 0, sticky = "e"
    )
    ##
    vari12 <- tclVar("")
    tkgrid(gme2d <- tkentry(textvariable = vari12,
                             parent = groupmeans, width = 20),
           row = 1, column = 1, sticky = "e", padx = 2
    )
    ##
    tkgrid(bgme2d <- tkbutton(text = gettext("help", domain = "R-MCPtests"),
                               parent = groupmeans, borderwidth = 0, image = tclvalue(imageinfo),
                              command = function(...) {
                                tkmessageBox(message = gettext("Enter the treatment levels. For example, for a character vector of four treatments: A, B, C, D. You do not need to use the concatenate function, i.e., c(). Nor will you need to use quotes between treatment levels.", domain = "R-MCPtests"))
                              }),
           row = 1, column = 2, sticky = "e"
    )
    ##
    ## DFerror
    tkgrid(la.dfe <- tklabel(text = gettext("DFerror:", domain = "R-MCPtests"),
                   parent = groupmeans),
           row = 2, column = 0, sticky = "e"
    )
    ##
    vari13 <- tclVar("")
    tkgrid(gme3d <- tkentry(textvariable = vari13,
                             parent = groupmeans, width = 20),
           row = 2, column = 1, sticky = "e", padx = 2
    )
    ##
    tkgrid(bgme3d <- tkbutton(text = gettext("help", domain = "R-MCPtests"),
                               parent = groupmeans, borderwidth = 0, image = tclvalue(imageinfo),
                              command = function(...) {
                                tkmessageBox(message = gettext("Enter the value of the degrees of freedom of the mean square error of the experiment model. The value is numeric.", domain = "R-MCPtests"))
                              }),
           row = 2, column = 2, sticky = "e"
    )
    ##
    ## MSerror
    tkgrid(la.mse <- tklabel(text = gettext("MSerror:", domain = "R-MCPtests"),
                   parent = groupmeans),
           row = 3, column = 0, sticky = "e"
    )
    ##
    vari14 <- tclVar("")
    tkgrid(gme4d <- tkentry(textvariable = vari14,
                             parent = groupmeans, width = 20),
           row = 3, column = 1, sticky = "e", padx = 2
    )
    ##
    tkgrid(bgme4d <- tkbutton(text = gettext("help", domain = "R-MCPtests"),
                               parent = groupmeans, borderwidth = 0, image = tclvalue(imageinfo),
                              command = function(...) {
                                tkmessageBox(message = gettext("Enter the value of the mean square error of the experiment model. The value is numeric.", domain = "R-MCPtests"))
                              }),
           row = 3, column = 2, sticky = "e"
    )
    ##

    ## Replication
    tkgrid(la.rep <- tklabel(text = gettext("Replication:", domain = "R-MCPtests"),
                   parent = groupmeans),
           row = 4, column = 0, sticky = "e"
    )
    ##
    vari15 <- tclVar("")
    tkgrid(gme5d <- tkentry(textvariable = vari15,
                             parent = groupmeans, width = 20),
           row = 4, column = 1, sticky = "e", padx = 2
    )
    ##
    tkgrid(bgme5d <- tkbutton(text = gettext("help", domain = "R-MCPtests"),
                               parent = groupmeans, borderwidth = 0, image = tclvalue(imageinfo),
                              command = function(...) {
                                tkmessageBox(message = gettext("Enter the value of number of replications of the treatments. The value is numeric and if the data is unbalanced use the harmonic mean of the replications. Inserted all the arguments above, click on the 'Calculate' button", domain = "R-MCPtests"))
                              }),
           row = 4, column = 2, sticky = "e"
    )
    ##
    ## Options 'Graphic Parameters'

    # Option 'Graphic Parameters'
    ###############################
    tkpack(frame.plot <- ttklabelframe(text = gettext("Graphic Parameters",
                                                      domain = "R-MCPtests"),
                                       parent = child21.group2, style = "Toolbar.TLabelframe"),
           expand = TRUE, fill = "both"
    )
    ## Color
    tkgrid(tklabel(text = gettext("Color:", domain = "R-MCPtests"),
                   parent = frame.plot),
           row = 0, column = 0, sticky = "e", pady = 3
    )
    ##
    vari16 <- tclVar("")
    tkgrid(ggraf1d <- tkentry(textvariable = vari16,
                               parent = frame.plot, width = 20),
           row = 0, column = 1, sticky = "e", padx = 2
    )
    ##
    tkgrid(bgraf1d <- tkbutton(text = gettext("help", domain = "R-MCPtests"),
                                parent = frame.plot, borderwidth = 0, image = tclvalue(imageinfo),
                               command = function(...) {
                                 tkmessageBox(message = gettext("Enter the color name of the chart bars. Place quotation marks around the color names. For example, if you want the color red for the bars, use 'red'. For functions, quotation marks are not necessary, for example gray.colors() function.", domain = "R-MCPtests"))
                               }),
           row = 0, column = 2, sticky = "e"
    )
    ##

    ## Horizontal
    tkgrid(tklabel(text = gettext("Horizontal:", domain = "R-MCPtests"),
                   parent = frame.plot),
           row = 1, column = 0, sticky = "e", pady = 3
    )
    ##
    vari17 <- tclVar("")
    tkgrid(ggraf2d <- tkentry(textvariable = vari17,
                               parent = frame.plot, width = 20),
           row = 1, column = 1, sticky = "e", padx = 2
    )
    ##
    tkgrid(bgraf2d <- tkbutton(text = gettext("help", domain = "R-MCPtests"),
                                parent = frame.plot, borderwidth = 0, image = tclvalue(imageinfo),
                               command = function(...) {
                                 tkmessageBox(message = gettext("Choose horizontal or vertical bars (FALSE or TRUE).", domain = "R-MCPtests"))
                               }),
           row = 1, column = 2, sticky = "e"
    )
    ##
    ## Axes
    tkgrid(tklabel(text = gettext("Axes:", domain = "R-MCPtests"),
                   parent = frame.plot),
           row = 2, column = 0, sticky = "e", pady = 3
    )
    ##
    vari18 <- tclVar("")
    tkgrid(ggraf3d <- tkentry(textvariable = vari18,
                               parent = frame.plot, width = 20),
           row = 2, column = 1, sticky = "e", padx = 2
    )
    ##
    tkgrid(bgraf3d <- tkbutton(text = gettext("help", domain = "R-MCPtests"),
                                parent = frame.plot, borderwidth = 0, image = tclvalue(imageinfo),
                               command = function(...) {
                                 tkmessageBox(message = gettext("Enter the axes. Separate them by semicolons. To add the names on the X and Y axes: xlab = 'Label X-axix'; ylab = 'Label Y-axis'.", domain = "R-MCPtests"))
                               }),
           row = 2, column = 2, sticky = "e"
    )
    ##
    ## Plot button
    updateplot <- function(...) {
      if (tclvalue(vari3) == gettext("Model", domain = "R-MCPtests")) {
        objtreat <- as.factor(envMCP$dat[,tclvalue(vari6)])
      }
      if (tclvalue(vari3) == gettext("Response variable", domain = "R-MCPtests")) {
        objtreat <- as.factor(envMCP$dat[,tclvalue(vari8)])
      }
      if (tclvalue(vari3) == gettext("Averages", domain = "R-MCPtests")) {
        trat <- strsplit(tclvalue(vari12), split = ",", perl = TRUE)[[1]]
        trat <- as.factor(trat)
        objtreat <- trat
      }
      xlab <- NULL; ylab <- NULL
      eval(parse(text = tclvalue(vari18)))
      hor <- if (tclvalue(vari17) == "") FALSE else eval(parse(text = tclvalue(vari17)))
      color <- if (tclvalue(vari16) == "") heat.colors(length(levels(objtreat))) else {
        eval(parse(text = tclvalue(vari16)))
        # if (any(c(grep("\'", tclvalue(vari16)), grep("\"", tclvalue(vari16))) == 1)) {
        #   if (length(grep("\'", tclvalue(vari16))) > 0) {
        #     gsub("\'", "", tclvalue(vari16))
        #   }
        #   if (length(grep("\"", tclvalue(vari16))) > 0) {
        #     gsub("\"", "", tclvalue(vari16))
        #   }
        # } else{
        #   eval(parse(text = tclvalue(vari16)))
        # }
      }
      #Update plots
      sapply(as.character(tkwinfo("children", frame.graf)),
             function(W) tcl("destroy", W))

      #Frame label
      tkpack(framelabel <- ttkframe(frame.graf, style = "Toolbar.TFrame",
                                    padding = c(3,3,3,3)),
             side = "top", fill = "x", expand = TRUE)

      # Frame label aux
      tkpack(framelabelaux <- ttkframe(framelabel, style = "Toolbar.TFrame",
                                       padding = c(3,3,3,3)),
             side = "top", fill = "x", expand = TRUE)

      #Frame plot
      tkpack(frameplot <- ttkframe(frame.graf, style = "Toolbar.TFrame",
                                   padding = c(3,3,3,3), height = 360),
             side = "top", fill = "both", expand = TRUE)

      # Creating an image in the tcltk
      fator <- 1
      #Plot
      ##
      plottopwinMCP <- tkplot(parent = frameplot,
                              function(...) {
                                MCPtests::MCPbarplot(envMCP$results, col = color, horiz = hor, xlab = xlab, ylab = ylab)
                              }, hscale = fator, vscale = fator)
      # Auxiliar function
      f <- function(...) {
        fatorup <- as.numeric(tclvalue("fator"))
        if (fator != fatorup) {
          fator <<- fatorup
          tkrreplot(plottopwinMCP, hscale = fator, vscale = fator)
        }
      }
      # Button save as...
      #graphics.off() # Erasing All Graphics Devices
      dispplot <- NULL
      bsaveas <- tkbutton(parent = framelabelaux, text = gettext("Save as...", domain = "R-MCPtests"),
                          borderwidth = 0, underline = 0,
                          image = tclvalue(imagesaveas), compound = "top",
                          command = function(...){
                            dispplot <<- grDevices::dev.new(noRStudioGD = TRUE) # New device plot
                            #grDevices::dev.new(noRStudioGD = TRUE) # New device plot
                            if (tclvalue(vari3) == gettext("Model", domain = "R-MCPtests")) {
                              objtreat <- as.factor(envMCP$dat[,tclvalue(vari6)])
                            }
                            if (tclvalue(vari3) == gettext("Response variable", domain = "R-MCPtests")) {
                              objtreat <- as.factor(envMCP$dat[,tclvalue(vari8)])
                            }

                            if (tclvalue(vari3) == gettext("Averages", domain = "R-MCPtests")) {
                              # Treatment levels
                              trat <- strsplit(tclvalue(vari12), split = ",", perl = TRUE)[[1]]
                              trat <- as.factor(trat)
                              objtreat <- trat
                            }
                            xlab <- NULL; ylab <- NULL
                            eval(parse(text = tclvalue(vari18)))
                            hor <- if (tclvalue(vari17) == "") FALSE else eval(parse(text = tclvalue(vari17)))
                            if (tclvalue(vari16) == "") {
                              color <- heat.colors(length(levels(objtreat)))
                            } else {
                              if (any(c(grep("\'", tclvalue(vari16)), grep("\"", tclvalue(vari16))) == 1)) {
                                if (length(grep("\'", tclvalue(vari16))) > 0) {
                                  color <- gsub("\'", "", tclvalue(vari16))
                                }
                                if (length(grep("\"", tclvalue(vari16))) > 0) {
                                  color <- gsub("\"", "", tclvalue(vari16))
                                }
                              } else{
                                color <- eval(parse(text = tclvalue(vari16)))
                              }
                            }
                            MCPtests::MCPbarplot(envMCP$results, col = color, horiz = hor, xlab = xlab, ylab = ylab)
                          })
      tkpack(bsaveas, side = "left")

      # Frame label aux2
      tkpack(framelabelaux2 <- ttkframe(framelabelaux, style = "Toolbar.TFrame",
                                        padding = c(3,3,3,3)),
             side = "top", fill = "x", expand = TRUE)
      #Label
      scalelabel <- tklabel(parent = framelabelaux2, text = gettext("Scale of Plot", domain = "R-MCPtests"))
      # Scale of plot
      s <- tkscale(framelabelaux2, command = f, from = 1, to = 3.00, variable = "fator",
                   showvalue = TRUE, resolution = 0.05, orient = "horiz", borderwidth = 0)

      # Drawing the widget
      tkpack(s)
      tkpack(scalelabel, before = s)
      tkpack.configure(s, fill = "both", expand = TRUE, side = "bottom", ipady = 5)
      #tkpack.configure(scalelabel, side = "bottom", anchor = "s")



      #tkgrid(plottopwinMCP,  row = 0, column = 0, sticky = "news")
      tkpack(plottopwinMCP,  fill = "both", expand = TRUE)
      addScrollbars(frameplot, plottopwinMCP)
    }
    tkpack(plot_button <- ttkbutton(text = gettext("Update plot",
                                                   domain = "R-MCPtests"),
                                    parent = child21.group2,
                                    command = updateplot),
           side = "bottom",
           anchor = "s",
           expand = TRUE,
           fill = "x"
    )
    tkbind(topwinMCP, "<Alt-Return>", updateplot)

    ############
    # Working Q4
    ############
    tkpack(frame.graf <- ttklabelframe(text = gettext("Plot",
                                                      domain = "R-MCPtests"),
                                       parent = child22.group2,
                                       style = "Toolbar.TLabelframe"
                                       ),
           expand = TRUE, fill = "both"
    )
    tkpack(tklabel(parent = frame.graf, image = "::image::MCP2"),
           expand = TRUE, fill = "both")
    tkpack(tklabel(parent = frame.graf, foreground = "blue",
                   text = gettext("Democratizing statistical knowledge to the world!", domain = "R-MCPtests")),
           side = "bottom", anchor = "center"
           )


    # Activate GUI
    finish <- tclServiceMode(oldmode)




    # Mensagem de protecao para nao fechar o programa por engano
    tkwm.protocol(topwinMCP, "WM_DELETE_WINDOW", function(){
      response <- tkmessageBox(icon = "question",
                               message = gettext("Really close?", domain = "R-MCPtests"),
                               type = "yesno" ,
                               parent = topwinMCP)
      tkdestroy(topwinMCP) # Caso contrario feche
    }
    )
    tkwm.deiconify(topwinMCP)
    #tkwm.state(topwinMCP, "zoomed")
  }
  if (gui == FALSE) {
    response <- tk_messageBox(
      title = gettext("Tell me something:", domain = "R-MCPtests"),
      message = gettext("Do you want to use the GUI for the package?", domain = "R-MCPtests"),
      icon = "question",
      type = "yesno"
    )
    if (response == "yes") {
      guiMCP(gui = TRUE)
    }
    if (response == "no") {
      tk_messageBox(message = gettext("Use the MRtest function! For help, use ?MRtest", domain = "R-MCPtests"))
      ?MCPtests::MCPtest
    }

  }

}
