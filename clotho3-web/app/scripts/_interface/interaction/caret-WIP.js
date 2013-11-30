/*****************
 text caret selection and movement

 - see framework in references folder

 inspiration:
 - https://github.com/akatov/angular-contenteditable/blob/master/angular-contenteditable.js
 - https://github.com/DrPheltRight/jquery-caret/blob/master/jquery.caret.js
 - http://plugins.jquery.com/caret/
 - https://code.google.com/p/rangyinputs/
 - http://stackoverflow.com/questions/1181700/set-cursor-position-on-contenteditable-div?rq=1
 ******************/

//note - jQuery dependence
angular.module('clotho.interface').service('$caret', function($log) {
	var functions = {};

	//handle browser idiosyncrasies, and element type (input, contenteditable, etc.)
	//todo - account for multiple nodes -- check in whole element (look for parent e.g. plasmid editor)

	function caretTo(el, index) {
		if (el.createTextRange) {
			var range = el.createTextRange();
			range.move("character", index);
			range.select();
		} else if (el.selectionStart != null) {
			el.focus();
			el.setSelectionRange(index, index);
		}
	}



	//adds DOM node to mark selection
	functions.savePos = function() {
		//todo
	};

	//removes DOM node and restores caret / selection to original position (of savePos)
	functions.restorePos = function() {
		//todo
	};



	//works
	//future - get cursor location, and reset to there rather than end
	//note - expects DOM element, not jQuery
	functions.setEndOfContenteditable = function(contentEditableElement)
	{
		var range,selection;
		if(document.createRange)//Firefox, Chrome, Opera, Safari, IE 9+
		{
			range = document.createRange();
			range.selectNodeContents(contentEditableElement);
			range.collapse(false);
			selection = window.getSelection();
			selection.removeAllRanges();
			selection.addRange(range);
		}
		else if(document.selection)//IE 8 and lower
		{
			range = document.body.createTextRange();
			range.moveToElementText(contentEditableElement);
			range.collapse(false);
			range.select();
		}
	};

	functions.getPos = function(el) {
		//todo - account for multiple nodes -- check in whole element (look for parent e.g. plasmid editor)

		var saved;

		if(window.getSelection)//HTML5
		{
			if (el.contenteditable == 'true') {
				var range = window.getSelection().getRangeAt(0),
					start = range.cloneRange(),
					end = range.cloneRange();
				start.collapse(true);
				end.collapse(false);
				saved = start; //todo
			} else {
				saved =  el.selectionStart;
			}
		}
		else if(document.selection)//IE<9
		{
			if (el.contenteditable == 'true') {
				var r1 = window.getSelection().getRangeAt(0),
					r2 = r1.cloneRange();
				r2.moveToElementText(el);
				r2.setEndPoint('EndToEnd', el);
				saved = r2.text.length;
			} else {
				saved = document.selection.createRange();
				//todo - check this one
			}
		} else {
			saved = 0;
		}

		$log.log(saved);

		return saved;

	};

	functions.setPos = function(el, index) {
		switch (index) {
			case angular.isEmpty(index) : {
				index = 0;
			}
			case 'end' : {
				index = el.val().length();
			}
			case 'start' : {
				index = 0;
			}
			default : {}
		}

		//check type - contenteditable or input
	};

	/**
	 * @return {array} [start, end]
	 */
	functions.getSel = function() {
		//todo - see plasmid module
	};


	functions.setSel = function(el, start, end) {
		//todo
	};


	functions.insert = function (node, content) {
		//todo - see plasmid module
	};

	return functions;
});