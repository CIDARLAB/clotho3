/**
 * @fileOverview  Searchbar script for Clotho project cl01
 * @author werner@bussedesign.com
 * @version 1
 * @requires jQuery 1.8.3
 * @reguires jQuery UI 1.10
 */

(function ($) {
    $(function () {
        // we are using a closure to limit variable scope to this function only.
        // without this we might have a problem with firebug showing the html that is loaded after an ajax call

        /**
         *  Render the search item details
         *  the function uses two data types and an associated json object
         *  @param {Object} jsonObj, the associated data object
         *  @param {String} dataLabel, the name of the object or the author
         *  @param {String} dataTyle, either content or author
         *  @returns {Array} array with members infoBlockContentHeader(String), listItem(String)
         */
        function renderSearchItemDetails(jsonObj, dataLabel, dataType) {

            var infoBlockContentHeader, listItem;
            
            if (dataType == 'content') {
                infoBlockContentHeader = "<div class='contentHeaderWrap'>";
                infoBlockContentHeader += "<img alt='' src='" + jsonObj.header.icon + "' />";
                infoBlockContentHeader += "<span class='headerText'>" + jsonObj.header.title + "</span>";
                infoBlockContentHeader += "<span class='rateDown'>" + jsonObj.header.rateDown + "</span>";
                infoBlockContentHeader += "<span class='rateUp'>" + jsonObj.header.rateUp + "</span>";
                infoBlockContentHeader += "</div>";
                listItem = "<ul class='contentList'>";
                for (var key in jsonObj.content[0]) {

                    listItem += "<li>";
                    listItem += "<span class='label'>" + key + ":</span>";
                    listItem += "<span class='value'>" + jsonObj.content[0][key] + "</span>";
                    listItem += "</li>";
                }
                listItem += "</ul>";
            } else if (dataType == 'author') {
                infoBlockContentHeader = "<div class='authorHeaderWrap'>";
                infoBlockContentHeader += "<img class='authorAvatar' alt='' src='" + jsonObj.avatar + "' />";
                infoBlockContentHeader += "<span class='headerText'>" + dataLabel + "</span>";
                infoBlockContentHeader += "<a class='authorEmail' href='mailto:" + jsonObj.email + "'>" + jsonObj.email + "</a>";
                infoBlockContentHeader += "</div>";
                listItem = "<p class='authorDesc'>" + jsonObj.description + "</p>";
            } // end of switch statement

            // return the rendered header content and the content area
            return [infoBlockContentHeader, listItem];
        } // end renderSearchItemDetails()





        /**
         *  Build the command string
         *  the function uses two data types and an associated json object
         *  @param {Object} jsonObj, the command object
         *  @returns {String} command string
         */
        function buildCommandString(jsonObj) {
            // build the command arguments
            var cmdArgs = '';
  
            // loop through the jason objects and construct the arguments
            for (var i = 0; i < jsonObj.commandArguments.length; i++) {
                cmdArgs += "'" + jsonObj.commandArguments[i].linkText + "'";
                if ((i + 1) < jsonObj.commandArguments.length) {
                    cmdArgs += ",";
                }
            }
            // return command string the command e.g. commandName(argument1, argument2,...)
            return jsonObj.commandName + "(" + cmdArgs + ")";
        } // end buildCommandString()





        /**
         * error message pane handling
         * @param {String} msg, the message to be displayed
         */
        function showMessage(msg) {
            msg = "<p>" + msg + "</p>";
            $('#msgPane').append(msg);
            $('#messages, #msgPane').addClass('active').fadeIn('slow');
        }

        function hideMessage() {          
            $('#messages, #msgPane').removeClass('active').fadeOut('slow', function(){
              $('#msgPane').find('p').remove();
            });        
        }
        
        
        $('#msgPane').find('.close').click(function () {
            hideMessage();
        });        
        
        
        
        
        /**
         *  submit search term
         *  this function may be invoced by hitting the ENTER key or by 
         *  clicking either the Submit or the EXECUTE button.
         *  it updates the activities log with the submitted search term and the 
         *  corresonding server response
         */
        function submitSearchTerm() {
        
            var searchSubmitURL = "http://www.bussedesignstage.com/ajax/formSubmitted.php";
            //var searchSubmitURL = "http://cl01.1/scripts/php/formSubmitted.php";
        
            // close the dropdown first
            $('ul.ui-autocomplete').fadeOut('slow');
            
            // validate that htere is something in the field
            if($("#searchBar").val() == '') {
              return
            }
            
            // update activity log with search term submission
            jsonObj = ({ 
                        msgBody: $("#searchBar").val(), // get the value of the search input field
                        msgType: 'regular', 
                        undoURL: ''
                      });
            updateActivityLog (jsonObj, 'client');
            
            // submit the search term via an ajax call to stay on page
            $.ajax({
                type: 'POST',
                url: searchSubmitURL,
                dataType:"jsonp",
                data: {  searchTerm: $SearchBar.val() },
                success: function(data){
                // update the activity log with server response
                updateActivityLog (data, 'server');
                }
            });
            
        } // end submitSearchTerm

        /**
         *  this is the data source for all ajax requests. data is deliverred in response to a url attribute
         *  e.g. http://www.bussedesignstage.com/ajax/remote.php?term=[searchterm]
         *    term = returns json object with search suggestions
         *    content = returns json object with block content
         *    author = return json object with author data
         *
         */
        var dataSource = "http://www.bussedesignstage.com/ajax/remote.php";
        
        /**
         * initialize searchBar for autocomplete
         */

        var $SearchBar =  $("#searchBar");
        var $loading = $("#loading");

        $SearchBar.autocomplete({
            source:function (request, response) {
                // callback to connect to the source
                $.ajax(
                    // the request object
                    {
                        url:dataSource + "?term=" + request.term,
                        dataType:"jsonp",
                        success:function (data) {
                            response(data);
                        }, // end success
                        timeout:10000, // 2 seconds timeout
                        beforeSend:function () {
                            $(".searchbarWrapper").append($loading);
                            $loading.show();
                            hideMessage();
                        },
                        complete:function () {
                            $loading.hide();
                        },
                        error:function () {
                            response({});
                            showMessage('The system is taking longer than expected to respond. Please try executing your command again in a few minutes.');
                        } // end error
                    })
            }, // end source
            delay:500,
            minLength:1

        })// end autocomplete






        /**
         * overwrite the single key/value pair default render method
         */
            .data("uiAutocomplete")._renderItem = function (ul, item) {

            //build the search suggestions drop-down list item from the inside out

            // build the link content
            // prepare the byline if there is one
            if (item.hasOwnProperty('byline')) {
                var linkContent = "<a class='searchSuggestion'>" + item.searchTerm + "<span class='byline'>" + item.byline + "</span></a>";     
            } else {
                linkContent = "<a class='searchSuggestion'>" + item.searchTerm + "</a>";
            }


            // construct a list item and insert the search suggestion
            // check for infoPane present. If yes use <li> with class hasInfoPane, that allows us to color the <li> differently
            // indicating that there is more
            if (item.hasOwnProperty("infoPane")) {
            
                var hasInfoPane = true;
                
                // buildListItem is used to assemble the list item including the info pane.
                // it is then appendedTo the <ul>
                var $buildListItem = $("<li class='hasInfoPane'>").append(linkContent);
                
            } else {
                $buildListItem = $("<li>").append(linkContent);
            }
            
            
            
            
            
            
            /**
             * construct the info pane if there is one
             */

            /**
             ---- .infoPane ----------------------------------------------------------------
             |                                                                             |
             |  ---- .header ------------------------------------------------------------  |
             |  |              Title, Command structure, Execute Button                 |  |
             |  -------------------------------------------------------------------------  |
             |                                                                             |
             |  ---- .content ----------------------------------------------------------   |
             |  |                                                                       |  |
             |  |  ---------- . blockHeader --------------------------------------      |  |
             |  |  |                                                             |--    |  |
             |  |  |      One block header per content block                     | |    |  |
             |  |  |      title button and author button                         | |--  |  |
             |  |  |                                                             | | |  |  |
             |  |  --------------------------------------------------------------- | |  |  |
             |  |    --------------------------------------------------------------- |  |  |
             |  |      ---------------------------------------------------------------  |  |
             |  |                                                                       |  |
             |  |  ------------ .blockContent ----------------------------------------  |  |
             |  |  |                                                                 |  |  |
             |  |  |  ---------- .infoBlockContentHeader --------------------------  |  |  |
             |  |  |  |                                                           |  |  |  |
             |  |  |  |  --- .contentHeaderWrap / .authorHeaderWrap ------------  |  |  |  |
             |  |  |  |  |                                                     |  |  |  |  |
             |  |  |  |  -------------------------------------------------------  |  |  |  |
             |  |  |  -------------------------------------------------------------  |  |  |
             |  |  |                                                                 |  |  |
             |  |  |  -------- .infoBlockContentBody ------------------------------  |  |  |
             |  |  |  |                                                           |  |  |  |
             |  |  |  |                                                           |  |  |  |
             |  |  |  -------------------------------------------------------------  |  |  |
             |  |  -------------------------------------------------------------------  |  |
             |  -------------------------------------------------------------------------  |
             -------------------------------------------------------------------------------

             */
            if (hasInfoPane) {

                // build the infoPane
                var infoPane = "<div class='infoPane'>";
                infoPane += "<div class='inner'>";
                
                
                
                
                // build the header
                var infoPaneHeader = "<div class='header'>";
                
                    // insert the title
                    infoPaneHeader += "<h2>" + item.searchTerm + "</h2>";
                    
                    
                    // if there is a command string to be displayed, build command string
                    if (item.infoPane.hasOwnProperty('command')) {
                        var commandString = buildCommandString(item.infoPane.command);
                        // construct the command e.g. commandName(argument1, argument2,...) and insert into header
                        infoPaneHeader += "<p class='command'>" + commandString + "</p>";                   
                    } // end build command string
    
    
                    // if there is an execute button to be displayed, build button
                    if (item.infoPane.hasOwnProperty('executeButton')) {
                        infoPaneHeader += "<a class='execute button' href='" + item.infoPane.executeButton.linkURL + "'>" + item.infoPane.executeButton.linkText + "</a>";
                    } // end build execute button
               
               
                infoPaneHeader += "</div>";  // close header
                
                
                // add the header to the info pane
                infoPane += infoPaneHeader;




                // create the info pane content
                if (item.infoPane.hasOwnProperty('infoContent')) {

                    var infoPaneContent = "<div class='content'>";
                    var infoBlockHeaders = "";
                    var i = 0;

                    // cycle through the infoBlock objects and build the headers
                    for (key in item.infoPane.infoContent) {
                    
                        // first header button to be active by default
                        if (i == 0) {
                            var state = 'active';
                        } else {
                            state = '';
                        }
                        
                        // build the header buttons
                        infoBlockHeaders += "<div class='blockHeader'>";
                        infoBlockHeaders += "<a class='infoBlockHeaderButton title " + state + "'>" + item.infoPane.infoContent[i].header.title + "</a>";
                        infoBlockHeaders += "<a class='infoBlockHeaderButton author'><img class='authorAvatar' src='" + item.infoPane.infoContent[i].header.author.avatar + "' alt='' /><span>Author:</span>" + item.infoPane.infoContent[i].header.author.linkText + "</a>";
                        infoBlockHeaders += "</div>";
                        i++;
                    }
                    
                    // Add all available block headers to the content pane
                    infoPaneContent += infoBlockHeaders;      


                    // build the content block
                    var infoBlockContent = "<div class='blockContent'>";
                    
                    // add the info block content header
                    infoBlockContent += "<div class='infoBlockContentHeader'>";
                                       
                    // add first block content to content pane, all other blocks will be loaded later
                    // this block contains a sub header and the info
                    // renderSearchItemDetails(jsonObj, dataLabel, dataType)
                    var renderedItems = renderSearchItemDetails(item.infoPane.infoContent[0].content, item.infoPane.infoContent[0].header.title, "content");                    
                    infoBlockContent += renderedItems[0];
                    infoBlockContent += "</div>"; // close infoBlockContentHeader
                        
                    // insert the content
                    // first wrap list in div container
                    infoBlockContent += "<div class='infoBlockContentBody'>";
                    infoBlockContent += renderedItems[1];
                    infoBlockContent += "</div>"; // close infoBlockContentBody
                    
                    // insert the info block content into the info pane content
                    infoBlockContent += "</div>"; // close blockContent 
                                        
                    // Add content block to the content pane
                    infoPaneContent += infoBlockContent;
                                      
                    // add content pane to the info pane
                    infoPane += infoPaneContent;
                    
                }
                               
                infoPane += "</div>"; // close inner div
                infoPane += "</div>"; // close infoPane
                                
                // add the completed info pane to the list item
                $buildListItem.append(infoPane);

            }


            // add complete list item to the list
            return $buildListItem.appendTo(ul);

            // at this point we have build the complete search term suggestion list with all associated info panes

        }; // end rendering method






        /**
         * set some auto complete options
         */
         
        // place selected search term into search field
        $SearchBar.on("autocompleteselect", function (event, ui) {
            if (ui.item.hasOwnProperty("infoPane")) {
                // put command string into searchbar
                var commandString = buildCommandString(ui.item.infoPane.command);
                $SearchBar.val(commandString);
            } else {
                // no infopane = no command string
                $SearchBar.val(ui.item.searchTerm);
            }
            return false;
        });


        // place search term into search field on list item focus
        $SearchBar.on("autocompletefocus", function (event, ui) {
            if (ui.item.hasOwnProperty("infoPane")) {
                // put command string into searchbar
                var commandString = buildCommandString(ui.item.infoPane.command);
                $SearchBar.val(commandString);
                
            } else {
                // no infopane = no command string
                $SearchBar.val(ui.item.searchTerm);
            }
            return false;
        });
        
       

        // change click behavior for inside infoPane
        $ui_autocomplete = $(".ui-autocomplete");
        $ui_autocomplete.on("click", ".infoPane", function () {
            return false;
        });
      
      
        // allow click on author email to open email client
        $ui_autocomplete.on("click", ".authorEmail", function (event) {
            event.stopPropagation()
            return;
        });
        
        
        
        
        
        
        /**
         * click behavior for all infoBlockHeaderButtons
         */
         
        $ui_autocomplete.on("click", ".infoBlockHeaderButton", function () {
            // check if the button clicked is active already
            if ($(this).hasClass("active")) {
                return false;
            }


            // reset previous active button and set active class on clicked button
            $(".infoBlockHeaderButton").removeClass("active");
            $(this).addClass("active");


            // get data for this button
            if ($(this).hasClass('title')) {
                var dataType = "content";
                var dataLabel = $(this).html();
            }
            if ($(this).hasClass('author')) {
                dataType = "author";
                // extract just the name
                dataLabel = $(this).text().substring(7);
            }
            
            
            // store clicked DOM object for later use
            var $thisButton = $(this);


            // check if we already have the data
            // we have already gotten the data before and stored them in this dom element with the data attribute
            if ($thisButton.data('jsonData')) {

                // remove the old content data                
                // find the associated blockContent
                var $blockContent = $thisButton.parents('.content').children('.blockContent');
                $blockContent.fadeOut('slow', function () {
                    var $buildData = $thisButton.data('jsonData');
                    // render block content
                    // renderSearchItemDetails(jsonObj, dataLabel, dataType)
                    var renderedItems = renderSearchItemDetails($buildData, dataLabel, dataType);
                    // insert the new content header
                    $blockContent.find('.infoBlockContentHeader').html(renderedItems[0]);
                    // insert new content
                    $blockContent.find('.infoBlockContentBody').html(renderedItems[1]);
                    // make the data visible
                    $blockContent.fadeIn('slow');
                });
                return false;
            }
            
            

            // data has not been fetched, get it
            $.ajax(
                {
                    url:dataSource + "?" + dataType + "=" + dataLabel,
                    dataType:"jsonp",
                    success:function (data) {
                   
                        // store data in DOM object
                        $thisButton.data('jsonData', data);


                        // remove the old content data
                        // find the associated blockContent
                        var $blockContent = $thisButton.parents('.content').children('.blockContent');
                        $blockContent.fadeOut('slow', function () {
                        
                            // render block content
                            // renderSearchItemDetails(jsonObj, dataLabel, dataType) 
                            var renderedItems = renderSearchItemDetails(data, dataLabel, dataType);
                            // insert the new content header
                            $blockContent.find('.infoBlockContentHeader').html(renderedItems[0]);
                            // insert new content
                            $blockContent.find('.infoBlockContentBody').html(renderedItems[1]);
                            // make the data visible
                            $blockContent.fadeIn('slow');
                            
                        });
                        
                    }, // end success
                    beforeSend:function () {
                        $thisButton.append($loading);
                        $loading.show();
                    },
                    complete:function () {
                        $("#loading").hide();
                    },
                    timeout:5000, // 5 seconds timeout
                    error:function () {
                    var $blockContent =  $("#blockContent");
                        $blockContent.fadeOut('slow', function () {
                            $blockContent.find('.infoBlockContentHeader').html("");
                            $blockContent.find('.infoBlockContentBody').html("<p class='timeout'>The system is taking longer than expected to respond. Please try again in a few minutes.'</p>");
                            $blockContent.fadeIn('slow');
                        });
                    } // end error
                });

            return false;
        });


        /**
         *  Searchbar submissions
         *  this is where we submit all searches
         *  Since we will stay on the page we intercept all regular submissions and send them via an ajax call
         *
         *  there are three ways to submit a search
         *  1 via the keyboard ENTER key
         *  2 via the submit button in the searchbar
         *  3 via the EXECUTE button in the command info pane
         *
         */
         
        // click behavior for execute button
        $ui_autocomplete.on("click", ".execute", function () {
            submitSearchTerm();
            return false;
        }); 
     
        // click behavior for submit button
        $('#searchBarSubmit').click(function () {
            submitSearchTerm();
            return false;
        });
        
        // submit via ENTER key
        $('html').keydown(function (event) {
            if(event.keyCode == 13) {
                submitSearchTerm();
                return false;
            }    
        });
        
        

    });
})(jQuery);