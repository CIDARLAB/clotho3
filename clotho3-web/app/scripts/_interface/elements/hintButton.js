angular.module('clotho.interface').directive('hintButton', function() {
	return {
		scope: {
			hint : '@hintButton'
		},
		replace: true,
		template : '<button class="btn" popover="{{ hint }}" popover-trigger="mouseenter" popover-placement="left"><i class="icon-info-sign"></i> Hint</button>',
		link: function(scope, element, attrs) {}
	}
});