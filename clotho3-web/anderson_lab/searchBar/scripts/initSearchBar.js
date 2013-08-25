/**
 * @fileOverview initSearchBar for  Clotho project cl01.1
 * @author werner@bussedesign.com
 * @version 1
 * @requires jQuery 1.8.3
 * @reguires jQuery UI 1.10
 */

(function ($) {
    $(function () {
    
        var thisFilename = "initSearchBar.js";
                  
        /**
         * get the attributes of the javasript file we have injected
         */            
        //var backgroundColor = typeof sbarAttr.backgroundColor !== 'undefined' ? sbarAttr.backgroundColor : 'transparent';
        var backgroundIsWhite = typeof sbarAttr.backgroundIsWhite !== 'undefined' ?  sbarAttr.backgroundIsWhite : 'true';
        if (backgroundIsWhite) {
            backgroundColor = "#ffffff";
        } else {
            backgroundColor = "#000000";
        }
        var showHelpPane = typeof sbarAttr.showHelpPane !== 'undefined' ?  sbarAttr.showHelpPane : 'true';
        var backgroundAlpha = typeof sbarAttr.backgroundAlpha !== 'undefined' ?  sbarAttr.backgroundAlpha : '1';
        var buttonAlpha = typeof sbarAttr.buttonAlpha !== 'undefined' ?  sbarAttr.buttonAlpha : '1';
        
        /**
         * insert the searchbar style sheet link tag into the head of the page
         */
        //$('head').append("<link href='http://cl01.1/searchBar/styles/normalize.css' rel='stylesheet'>");
        //$('head').append("<link href='http://cl01.1/searchBar/styles/searchbar.css' rel='stylesheet'>");
        //$('head').append("<link href='http://bussedesignstage.com/ajax/cl01.1/searchBar/styles/normalize.css' rel='stylesheet'>");
        //$('head').append("<link href='http://bussedesignstage.com/ajax/cl01.1/searchBar/styles/searchBar.css' rel='stylesheet'>"); 
         
      
        
        /**
         *  build the searchbar html
         */
        
        var searchbarHTML  = "<div id='searchBarBackground'></div><div class='ui-widget'>";
            searchbarHTML += "<form id='searchBarForm' method='post'>";
            searchbarHTML += "<div class='searchbarWrapper ui-front'>";
            searchbarHTML += "<input id='searchBar' name='searchBar' type='text' placeholder='Enter Command or Search Term' />";
            searchbarHTML += "<input id='searchBarSubmit' type='submit' value='submit' />";
            // add the loading indicator
            searchbarHTML += "<div id='loading'></div>";
            
            // this really should be injected later with all other help panes
            searchbarHTML += "<div class='helpPane'>";
            searchbarHTML += "<span class='close'>x</span>";
            searchbarHTML += "<div class='helpContent'>";
            searchbarHTML += "<p class='helpMessage'>You can input commands to execute in your work area, or enter a search to find existing lorem ipsum dolores.</p>";
            searchbarHTML += "</div>";
            searchbarHTML += "<div class='helpControl'>";
            searchbarHTML += "<label class='neverAgain'><input type='checkbox'>Don’t show me tips again.</label>";
            searchbarHTML += "</div>";
            searchbarHTML += "</div>";
            searchbarHTML += "</div>"; // end searchBarWrapper

            // activity log 
            searchbarHTML += "<div id='log'>";
            // the notification bubble
            searchbarHTML += "<div id='notify'></div>";
            searchbarHTML += "<a>log</a>";
            searchbarHTML += "<div id='logPane'><h2>Activity Log</h2><a class='close'>x</a>";
            searchbarHTML += "<div class='inner'><ul id='logContent'>";
            
            // build the time stamp   
            var msgDate = new Date();
            var timestamp = msgDate.getDate() + "/" + (msgDate.getMonth()+1) + "/" + msgDate.getFullYear() + "  " + msgDate.getHours() + ":" + msgDate.getMinutes() + ":" + msgDate.getSeconds();
            // initial session message
            searchbarHTML += "<li class='client'><p class='message'>Session started</p><p class='timestamp'>" + timestamp + "</p></li>";
            
            searchbarHTML += "</ul></div></div>";
            searchbarHTML += "</div>"; // end activity log

            searchbarHTML += "</form>";  // end searchBarForm
            searchbarHTML += "</div>";  // end ui-widget
          
            //communication links
            searchbarHTML += "<ul id='communications'>";
            searchbarHTML += "<li id='newPage'><a href='#' title='New Page'>&nbsp;</a></li>";
            searchbarHTML += "<li id='newWorkSpace'><a href='#' title='New Work Space'>&nbsp;</a></li>";
            searchbarHTML += "<li id='info'><a>&nbsp;</a>";
            searchbarHTML += "<ul>";
            searchbarHTML += "<li><a href=''>Show me how</a></li>";
            searchbarHTML += "<li><a href=''>About Clotho</a></li>";
            searchbarHTML += "<li id='tooltips'><a data-tooltip='off' href=''>Turn tooltips on</a><a data-tooltip='on' href=''>Turn tooltips off</a></li>";
            searchbarHTML += "</ul>";
            searchbarHTML += "</li>";
            searchbarHTML += "</ul>";
            
            
            /**
             * add the searchbar to the DOM
             */
            $searchBarContainer = $('#searchBarContainer');
            
            
            // we append the html with a callback after the div has been hidden to avoid a flash of raw html            
            $searchBarContainer.hide('fast', function() { 
                $(this).append(searchbarHTML)
            });
            
            // set the background color
            
            if(!(backgroundColor == 'none')) {
            
                // convert the hex value into rgba(r,g,b,alpha)
                // get the r value
                var r = parseInt(backgroundColor.slice(1, 3), 16);
                var g = parseInt(backgroundColor.slice(3, 5), 16);
                var b = parseInt(backgroundColor.slice(5), 16);
                // assemble the css string and add alpha value
                backgroundColor = "rgba(" + r + "," + g + "," + b + "," + backgroundAlpha + ");";
                // set background color
            
            }
            $searchBarContainer.css('background-color', backgroundColor); 
                        
            
            /**
             * add the error message pane to the DOM
             */
            var errorPaneHTML  = "<div id='msgPane'>";
                errorPaneHTML += "<span class='alert'></span>";
                errorPaneHTML += "<span class='close'></span>";
                errorPaneHTML += "</div>";
                
            $(errorPaneHTML).insertAfter($searchBarContainer);
            
            /**
             * now fade in the searchbar
             */    
            $searchBarContainer.fadeIn('slow', function() {
                    
                
                // set final opacity
                //$searchBarContainer.css('opacity', opacity);
                              
                
                /**
                 * search bar specific help pane
                 */
                 
                if(showHelpPane) { 
                    // debug only -- comment this out before deployment
                    localStorage.removeItem('hideTips');
                     
                     
                    // if we set the hide help panes before we will not show the panes 
                    if(!localStorage.getItem('hideTips')) { 
                        $('#searchBarContainer').find('.helpPane').fadeIn('slow');
                        var activeTimerSearchBar = setTimeout(function(){
                            $('#searchBarContainer').find('.helpPane').fadeOut('slow'); 
                        }, 5000);
                    }
                    
                    
                    // we are showing help panes, each help pane has a close icon to hide it
                    $('#searchBarContainer').find('.helpPane').find('.close').click(function() {
                        if(activeTimerSearchBar) {
                            clearTimeout(activeTimerSearchBar);
                        }
                        $(this).parents('.helpPane').fadeOut('slow');
                    });
                    
                    
                    // we offer user the choice to never show tips again.
                    // This choice is made persistent by using a flag in local storage
                    $('#searchBarContainer').find('.helpPane').find('.neverAgain input').click(function() {
                        $(this).parents('.helpPane').fadeOut('slow');
                        // Put a flag into storage
                        localStorage.setItem('hideTips', 'true');
                    });
                    
                    
                    // if user clicks in search field we hide the help pane
                    $('#searchBar').focus(function(){
                        if(activeTimerSearchBar) {
                            clearTimeout(activeTimerSearchBar);
                        }
                        $('#searchBarContainer').find('.helpPane').fadeOut('slow');
                    });
                  }
                  
                  /**
                   * toggle the tooltip menu item
                   */ 
                  $('#tooltips').find("a[data-tooltip='off']").hide();
                   
                  $('#tooltips').find('a').click(function() {
                    if($(this).data('tooltip') == 'on') {
                      $(this).hide();
                      $("#tooltips a[data-tooltip='off']").show();
                    }
                    if($(this).data('tooltip') == 'off') {
                      $(this).hide();
                      $("#tooltips a[data-tooltip='on']").show();
                    }
                    return false;
                  });
                  
                  
                  // attach click handler to the log button
                  $('#log > a').click(function(){
                      if($(this).parent().hasClass('active')) {
                          $(this).next().fadeOut('slow', function(){
                              // remove the active class from the log div
                              $(this).parents('#log').removeClass('active');
                          });  
                      
                      } else {
                          $(this).parent().addClass('active');
                          $(this).next().fadeIn('slow');
                          
                          // the latest is always on the bottom
                          $scrollDiv = $(this).parent().find('.inner');
                          $scrollDiv.scrollTop($scrollDiv[0].scrollHeight);
                          
                          // hide notification bubble
                          var $notify = $('#notify');
                          $notify.fadeOut('slow', function (){ 
                              $notify.html('').removeClass('error');
                          }); 
                       } 
                       // return false so we don't trigger the document click handler to close #log
                       return false;  
                  });
                  
                  // attach click handler to close icon for log pane
                  $('.close').click(function(){
                      // remove the log pane
                      $(this).parent().fadeOut('slow', function(){
                          // remove the active class from the log div
                          $(this).parents('#log').removeClass('active');
                      });
                  });
                  
                  // if click outside the activities log pane we'll close the log
                  $(document).click(function(){
                    if($('#log').hasClass('active')) {
                        $('#logPane').fadeOut('slow', function(){
                            // remove the active class from the log div
                            $('#log').removeClass('active');
                        });
                    }
                  });
                  
                  // but not if we click inside the log page
                  $('#logPane').click(function(){
                      // return false so we don't trigger the document click handler to close #log
                      return false; 
                  });
                  
                  
                  
                  
                  // load searchbar script AFTER searchbarHTML has been loaded and is visible
                  $.getScript("http://bussedesignstage.com/ajax/cl01.1/searchBar/scripts/searchBar.js");
                  //$.getScript("searchBar/scripts/searchBar.js");
                  
                  // load the activity log script
                  $.getScript("http://bussedesignstage.com/ajax/cl01.1/searchBar/scripts/activityLog.js");
                  //$.getScript("searchBar/scripts/activityLog.js");
                 
                   
              }); // end of fadeIn callback   
                
    });
})(jQuery);