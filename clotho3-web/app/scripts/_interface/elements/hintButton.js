angular.module('clotho.interface').directive('hintButton', function() {
	return {
		replace: true,
		template : '<button class="btn" popover="{{ hint }}" popover-trigger="mouseenter" popover-placement="left"><i class="glyphicon glyphicon-info-sign"></i> Hint</button>',
		link: function(scope, element, attrs) {
			scope.hint = attrs.hintButton;
		}
	}
});