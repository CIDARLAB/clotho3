$(document).ready(function(){
	$(".widget").each(function() {
		widget_space.activateWidget(this);
	});
});

var widget_space = new Object();

widget_space.show = function (json) {
	//TODO - make options for draggable, resizable, close-button, etc. ----- deal with postitioning and width/height as well
	
	//NOTES
	//not dealing with more right now
	//do we want better positioning logic? i.e. take size of widget into account
	var posx = typeof json.posx !== 'undefined' ? json.posx : 20;
	var posy = typeof json.posy !== 'undefined' ? json.posy : 20;
	var objWidth = typeof json.objWidth !== 'undefined' ? json.objWidth : 200;
	var objHeight = typeof json.objHeight !== 'undefined' ? json.objHeight : 230;
	
	//assign as element.style so that can be changed by jquery when interacted with (e.g. drag, resize)
	//default z-index to 1000 so shows up on top. will be reset to appropriate lower number when dragged
	var widget = $(
		'<div widget_id="' + json.widget_id + '" class="widget resizable" style="z-index: 1000; position: absolute; top: '+posy+'; left: '+posx+'; width: '+objWidth+'px; height: '+objHeight+'px;"> \n \
			<div class="widget_close_button"></div> \n \
			<div class="widget-container"> \n \
				' + json.content + ' \n \
			</div> \n \
		</div>');
	
	$('div[widget_id=' + json.parent_widget_id + ']').append(widget);
	eval(json.on_show);
	
	widget_space.activateWidget(widget);
}

//JCA:  major rewrite on 9/29/12
widget_space.update = function (updateObject) {
    try {
        //Pull out the bits so this isn't so cryptic
        var script = updateObject.script;
        var args = [];
        var tokens = '';
        
        //Bundle up all the arguments
        var dataLinks = updateObject.inputArgs; //List<token, uuid> for input data
        for(token in dataLinks) {
            var uuid = dataLinks[token];  //The uuid for a sharable in the Collector
            var sharable = collector[uuid];
            args.push(sharable);
            tokens += token;
//            tokens += ',';
        }
        //if tokens ends in , trim in (to be added, then unsilence above)
 
        //If an updater function is undefined, go ahead and define it
        if(updater===undefined) {
            var updater = null;
        }
        
        //Wrap up the update method in a function, then evaluate with the arguments
        var wrapper = "updater = function(" + tokens + ") {" + script + "};"
        eval(wrapper);
        updater(args);
    } catch(err) {
        alert(err);
    }


}

widget_space.activateWidget = function (activateMe) {
	//make draggable, resizable, and activate close button
	//later may want to check object to see if resizable etc and add those attributes only if desired
	$(activateMe).draggable({containment: "#widget_workspace", snap: true, stack: ".widget"}).resizable({containment: "#widget_workspace", alsoResize: $(activateMe).children(".widget-container")});
	$(activateMe).children(".widget_close_button").click(function() {
		$(this).parent().remove();
		libsend.call("unlinkWidget", $(this).parent().attr("widget_id"));
	});
}

/* probably want to make the close button a div that is inserted on
 * mouseover so that its not so intrusive - need to deal with the
 * overflow:hidden parameter of widgets
 */
