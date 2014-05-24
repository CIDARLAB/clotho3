angular.module('clotho.tokenizer')
/**
 * auto-grow directive by the "shadow" tag concept
 */
.directive('autoGrow', function($timeout) {
	return {
		link: function(scope, element, attr){
			var paddingLeft = element.css('paddingLeft'),
				paddingRight = element.css('paddingRight');

			var minWidth = 100;

			var $shadow = angular.element('<span></span>').css({
				'position': 'absolute',
				'top': '-10000px',
				'left': '-10000px',
				'fontSize': element.css('fontSize'),
				'fontFamily': element.css('fontFamily'),
				'white-space': 'pre'
			});
			element.after($shadow);

			var update = function() {
				var val = element.val()
						.replace(/</g, '&lt;')
						.replace(/>/g, '&gt;')
						.replace(/&/g, '&amp;')
					;

				// If empty calculate by placeholder
				if (val !== "") {
					$shadow.html(val);
				} else {
					$shadow.html(element[0].placeholder);
				}

				//extra to handle new letter before next $digest
				var calcWidth = $shadow[0].offsetWidth + 26;
				var newWidth = Math.max(calcWidth, minWidth) + "px";
				element.css('width', newWidth);
			};

			element.bind('keyup keydown blur', update);

			// Update on the first link
			// $timeout is needed because the value of element is updated only after the $digest cycle
			// TODO: Maybe on compile time if we call update we won't need $timeout
			$timeout(update);
		}
	}
});