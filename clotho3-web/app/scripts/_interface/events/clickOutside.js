// note - use click-outside-active attr.. should set that attr to false when take action with this directive
// note - requires jQuery has()
// future - to make more universal, see:
// https://raw.github.com/cowboy/jquery-outside-events/v1.1/jquery.ba-outside-events.js
// angular way of creating ng-directives
angular.module('clotho.interface').directive('clickOutside', function($document, $parse) {
	return function(scope, element, attr) {

		//todo - add watch for action
		var clickAction = $parse(attr['clickOutside']),
			active;

		scope.$watch(function () {
			return $parse(attr.clickOutsideActive)(scope);
		}, function (newval, oldval) {
			active = !!newval;
		});

		var handler = function (event) {

			if (active) {
				if (element.has(event.target).length == 0) {
					console.log('click captured');

					event.preventDefault();
					event.stopPropagation();
					scope.$apply( clickAction(scope, {$event:event}) );
				}
			}
		};

		$document.bind('click', handler);

		scope.$on('$destroy', function() {
			$document.unbind('click', handler);
		})
	}
});