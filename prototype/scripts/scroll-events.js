(function($)
{
    var latency = 200;
     
    $.fn.scrollstop = function(func, settings) {
		 
		 return this.each(function()
        {
            var elem = $(this);
 
            var movement = false;
             
            var displayTimer = null
            var movementTimer = null;
         
            //only need this piece if there is a time limit on when the mouse stop can occur.
            if(settings.latency != null)
            {
                var x = null;
                var y = null;
         
                var counter = 0;
                var counterMax = Math.ceil(settings.latency / 100);
                 
                elem
                .mouseover(function(e)
                {
                    movement = true;
                     
                    //check if movement has stopped to a maximum time of 100*counterMax, after that event will not run at all unless you re-mouseover
                    displayTimer = setInterval(function()
                    {
                        counter++;
                         
                        if(counter < counterMax)
                        {
                            if(!movement)
                            {
                                clearTimeout(displayTimer);//clear the timeout to avoid any funkiness
                                 
                                //set the coordinates for the event to the ones from the document event
                                e.pageX = x;
                                e.pageY = y;
                                 
                                func(e);
                            }
                            //else do nothing, just iterate
                        }else movement = false;//we can turn this off to avoid using the timeout in the mousemove
                    }, 100)
                })
            }
             
            elem
            .mouseout(function(e)
            {
                //kill this timers incase it's still running
                clearTimeout(displayTimer);
                clearTimeout(movementTimer);
                 
                counter = 0;//reset counter for when we mouseover again
                movement = false;//set movement back to false
                 
                settings.onMouseout(e);//call our mouseout
            })
            .mousemove(function(e)
            {
                x = e.pageX;
                y = e.pageY;
                 
                if(movement)//if we have moused over this will be on
                {
                    //clear timer and set again, this will determine our "stop" which will occur if mouse is in same position for the delayToStop time or more milliseconds
                    clearTimeout(movementTimer);
                    movementTimer = setTimeout(function()
                    {
                        movement = false;
                        if(settings.latency == null) func(e);
                    }, settings.delayToStop);
                }
                else
                {
                    settings.onStopMove(e);//call our mousemove - this is after the stop
                    movement = true;
                }
            });
        });
	}
	
})(jQuery);