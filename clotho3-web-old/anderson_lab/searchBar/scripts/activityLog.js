/**
 * @fileOverview Activity log script for Clotho project cl01
 * @author werner@bussedesign.com
 * @version 1
 * @requires jQuery 1.8.3
 * @reguires jQuery UI 1.10
 */


(function ($) {
    $(function () {
        // we are using a closure to limit variable scope to this function only.
        // without this we might have a problem with firebug showing the html that is loaded after an ajax call

        //var dataSource = "http://cl01.1/scripts/php/status.php";
        var dataSource = "http://bussedesignstage.com/ajax/status.php";
        
        // long polling  loop
        (function poll() {
        
            $.ajax({
                url:dataSource,
                dataType:"jsonp",
                success:function (data) {
                    // here is where we add the message to the activity log
                    updateActivityLog(data, 'server');
                }, // end success
                complete: poll
            });
                   
        })();
        
    });
})(jQuery);






/**
 *  Add a message to the activity log
 *  this GLOBAL function uses a json object and a location variable
 *  @param {Object} jsonObj, the data object
 *  @param {scope} location, either "server" or "client"
 *  @returns nothing
 */

var timerID;
var timerCount;
var notificationCycle = false;
var $notify = $('#notify');
var notificationsStack = [];


function updateActivityLog (jsonObj, location) {

  var activity  = "<li class='" +  location  + "'>";
  
  // check for msg type
  if(jsonObj.msgType == 'error') {
      activity += "<p class='message " + jsonObj.msgType + "'>" + "<span>Error:</span>" + jsonObj.msgBody + "</p>";
  } else {
      if(jsonObj.undoURL) {
          activity += "<p class='message " + jsonObj.msgType + "'>" + jsonObj.msgBody + "<span><a class='undo' href='" + jsonObj.undoURL + "'>(undo)</a></p>";
      } else {
          activity += "<p class='message " + jsonObj.msgType + "'>" + jsonObj.msgBody + "</p>";
      }
  }
  
      
  // build the time stamp   
  var msgDate = new Date();
  var timestamp = msgDate.getDate() + "/" + (msgDate.getMonth()+1) + "/" + msgDate.getFullYear() + "  " + msgDate.getHours() + ":" + msgDate.getMinutes() + ":" + msgDate.getSeconds();
      
      activity += "<p class='timestamp'>" + timestamp + "</p>";
      activity += "</li>";
      
          
  $('#logContent').append(activity);    



  // scroll the activities list to the bottom to show the latest update
  $scrollDiv = $('#logContent').parent();
  $scrollDiv.scrollTop($scrollDiv[0].scrollHeight);
  
  
  
  
  // if the activity log is open we just exit. If it is closed we flash the latest message in the notify bubble
  if(($('#log').hasClass('active'))) {
  
      // flash the notification stack and exit
      if (notificationsStack.length != 0) {
          notificationsStack = [];
      } 
      return;
      
    } else {
        // activity log is closed
        
        // return if client side message
        if(location === 'client') return;
        
        // push server message to notifications stack
        notificationsStack.push(jsonObj);
        
        // is timer running?
        // the timer will keep running until a notification cycle is completed. If a cycle is cancelled,
        // the timer will keep running and the timer count will be maneged by manageNotifications()
        if(timerID) return;
        
        // timer is not running, start a new notification cycle
        timerCount = 6;
        timerID = setInterval( 'manageNotifications()', 1000 ); // interval timing 1 sec       
       
  }
      
      
      
};  // end updateActivityLog()


/**
 *  manage the server notification
 *  a server message is shown during a notification cycle. If a new server message is received while
 *  a previous one is still showing, the previous one is canceled and the new is shown. 
 *  
 *  @returns nothing
 */
function manageNotifications() {
  
  // is there a notification on the stack?
  // if notification stack is empty then we just count down and finish a regular cycle
  if (notificationsStack.length == 0) {
  
      // decrement timer
      timerCount--;
      
      // if timer count is not 0, continue notification cycle
      if(timerCount != 0) return;
      
      // timer count is 0, end notification cycle
      
      // fadeout current message and remove message text
      $notify.fadeOut('slow', function (){ 
          $notify.html('').removeClass('error')
      });
      
      // clear timer
      clearInterval(timerID);
      timerID = '';
      
      // reset notificationCycle flag
      notificationCycle = false;
      
      return;
      
   } else {
       // yes there is a notification on stack
       
       // if notificationCycle is true then we are showing a message
       if(notificationCycle == true) {
       
           // end the cycle prematurely and start another one with the latest message
       
           // fadeout current message and remove message text
           $notify.fadeOut('slow', function (){ 
               $notify.html('').removeClass('error')
           });
           
           // reset notificationCycle flag
           notificationCycle = false;
           
                      
       } else {
           // start notification cycle
           // we got here with a notification on the stack but no active notification cycle
           
           // get server message from notifications stack
           var tempObj = notificationsStack[0];
           
           // remove server message from notifications stack
           notificationsStack = notificationsStack.slice( 1, notificationsStack.length );
           
           //Prepare the message
           if(tempObj.msgType == 'error') {
               $notify.addClass('error');
               var bubbleMsg = "Error: " + tempObj.msgBody;
           } else {
               var bubbleMsg = tempObj.msgBody;
           }
           
           // show server message
           $notify.append(bubbleMsg).fadeIn();
           
           // set notificationCycle flag
           notificationCycle = true;
           
           // start a new notification cycle
           timerCount = 6;
           
       }  
   }   
} // end manageNotifications











