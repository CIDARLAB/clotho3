$(document).ready(function() {
	$(document.body).append('<div id="tips_wrapper"></div>');
	$(document.body).append('<div id="lightbox"></div>');

	$("#theQuery-help").click(function(e) {
		e.preventDefault();
		help.display_all_tips();
	});
});

var help = new Object();

//later this will be sent from the server, or call help.add_tip_json
help.tips_json = JSON.parse('{"tips": [ {"selector" : "[name=Description]", "content" : "This is a simple tooltip to describe Description", "type" : "layover"} ] }');
help.rendered_json = JSON.parse('{}');
help.saved_tips = []; 
help.showingTips = false;

//old way of adding tip, easy constructor, adds JSON object to our array of all tips
help.add_tip = function(selector, content, type, pos_y, pos_x, height, width, additionalStyles) {
	var type = typeof type !== 'undefined' ? type : "tooltip",
	pos_x = typeof pos_x !== 'undefined' ? ', "pos_x" : "' +pos_x +'"': "";
	pos_y = typeof pos_y !== 'undefined' ? ', "pos_y" : "' +pos_y +'"' : "";
	height = typeof height !== 'undefined' ? ', "height" : "' +height +'"': "";
	width = typeof width !== 'undefined' ? ', "width" : "' +width +'"': "",
	additionalStyles = typeof additionalStyles !== 'undefined' ? ', "additionalStyles" : "' +additionalStyles +'"': "";
	help.tips_json.tips.push(JSON.parse('{"selector" : "' + selector + '", "content" : "' + content + '", "type" : "' + type + '"' + pos_x + pos_y + height + width + additionalStyles + '}'));
};

//pass a parsed json object to only include certain fields
help.add_tip_json = function(json_obj) {
	help.tips_json.tips.push(json_obj);
};

//defaults to tooltip if "type" not specified
help.generate_tip = function(json) {
	if (json.type == "layover") {
		return help.generate_layover(json);
	}
	else  {
		return help.generate_tooltip(json);
	}
}

//accepts a parsed JSON object
help.generate_tooltip = function(json) {
	var	selector = typeof json.selector !== 'undefined' ? json.selector : '',
	pos_x = typeof json.pos_x !== 'undefined' ? json.pos_x : $(selector).offset().left,
	pos_y = typeof json.pos_y !== 'undefined' ? json.pos_y : $(selector).offset().top + $(selector)[0].offsetHeight,
	height = typeof json.height !== 'undefined' ? ', "height: ' + json.height +'px;': "",
	width = typeof json.width !== 'undefined' ? ', "width: ' + json.width +'px;': "";
	
	var div = '<div class="tooltip" selector="'+escape(selector)+'" style="top: '+pos_y+'px; left: '+pos_x+'px"' +width + height+ json.additionalStyles + '>' + json.content + '<span class="arrow yellow topArrow arrowAtLeft"></span></div>';
	return div;
}

//accepts a parsed JSON object
help.generate_layover = function(json) {
	var selector = typeof json.selector !== 'undefined' ? json.selector : '',
	pos_x = typeof json.pos_x !== 'undefined' ? json.pos_x : $(selector).offset().left,
	pos_y = typeof json.pos_y !== 'undefined' ? json.pos_y : $(selector).offset().top,
	width = typeof json.width !== "undefined" ? json.width : $(selector)[0].offsetWidth,
	height = typeof json.height !== "undefined" ? json.height : $(selector)[0].offsetHeight;
	
	//horizontal align to middle if width is less than width of parent and pos_x not defined
	if (width != $(selector)[0].offsetWidth && typeof json.pos_x === 'undefined') { pos_x += (($(selector)[0].offsetWidth - width) / 2); }
	
	var div = '<div class="layovertip" selector="'+escape(selector)+'" style="top: '+pos_y+'px; left: '+pos_x+'px; width: '+width+'px; height: '+height+'px;' + json.additionalStyles + '"><div class="content">' + json.content + '</div></div>';
	return div;
}

//need a good way to access array by selector for this to work - UUIDs later should work...
help.display_tip = function(selector) {}

//need a good way to access array by selector for this to work - UUIDs later should work...
help.hide_tip = function(selector) {}

help.display_all_tips = function() {
	if (help.saved_tips != "") {
		help.saved_tips.fadeIn("fast");
	}
	else {
		for (i = help.tips_json.tips.length - 1; i >= 0; i -= 1) {
			var curTip = help.tips_json.tips[i];
			var newTip = $(help.generate_tip(curTip));	
			$("#tips_wrapper").append(newTip);
						
			//we should be able to do something like this - need to escape CSS selector attr properly
			//console.log(curTip.selector + " " + escape(curTip.selector));
			//$(".tooltip[selector="+escape(curTip.selector)+"]").animate({"backgroundColor" : "#FF9999"});
			
			//REFERENCE - help.rendered_json[escape(selector)] = [pos_x, pos_y, $(newTip)[0].scrollWidth, $(newTip)[0].scrollHeight];
			
			/*
			//TO DO - positioning to prevent overlap
			//check all previous tips' locations
			var collision = false;
			for (k = 0; k < i; k += 1) {
				var curSel = escape(help.tips_json.tips[k].selector);
				//this doesn't seem to work...
				var curElem = help.rendered_json[curSel];
				
				console.log(help.rendered_json);
				console.log(curSel);
				console.log(curElem);
				
				
			}
			
			//once in its spot, add it to the rendered list 
			help.rendered_json[escape(selector)] = [pos_x, pos_y, $(newTip)[0].scrollWidth, $(newTip)[0].scrollHeight];
			
			//console.log(help.rendered_json);
			*/
		}
	}
	help.showingTips = true;
	help.show_shade();
};

help.hide_all_tips = function() {
	help.saved_tips = $(".tooltip, .layovertip").fadeOut("fast");
	help.hide_shade();
	help.showingTips = false;
};

help.toggle_all_tips = function() {
	if (help.showingTips) {
		help.hide_all_tips();
	} else {
		help.display_all_tips();
	}
}

help.show_shade = function() {
	$(".shade-backdrop").remove();
	$("#page_wrapper").prepend('<div class="shade-backdrop"></div>');
	// clicking on anything that isn't a tip will hide the tips
	$(".shade-backdrop").one("mousedown", function(e) {
		e.preventDefault();
		help.hide_all_tips();
		help.hide_lightbox();
	});
}

help.hide_shade = function() {
	$(".shade-backdrop").fadeOut("fast");
}

//force tip re-rendering on window resize
//may want to change this to account for in-page changes too
$(window).resize(function() { 
	help.saved_tips = [];
});

/**************************************************
Lightbox
**************************************************/

help.show_lightbox = function(content) {
	//later, save content assuming we use this for more complex things
	var $lbox = $("#lightbox").empty()
	$lbox.html(content).css({"top": $(window).height()/2 -  $lbox.height()*0.5, "left": $(window).width()/2 -  $lbox.width()*0.5}).fadeIn();
	help.show_shade();
}

help.hide_lightbox = function () {
	$("#lightbox").fadeOut();
}