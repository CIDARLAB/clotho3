// todo - add activate attr, deregister listener based on it
// note - use click-outside-active attr
// note - requires jQuery has()
// future - to make more universal, see:
// https://raw.github.com/cowboy/jquery-outside-events/v1.1/jquery.ba-outside-events.js
// angular way of creating ng-directives
angular.module('clotho.interface').directive('clickOutside', function($document, $parse) {
	return function(scope, element, attr) {

		//todo - add watch for action
		var clickAction = $parse(attr['clickOutside']),
			active;

		attr.$observe('clickOutsideActive', function(value) {
			active = !!value;
		});

		var handler = function (event) {

			if (active) {
				event.preventDefault();
				event.stopPropagation();

				if (element.has(event.target).length == 0)
					scope.$apply( clickAction(scope, {$event:event}) );

				//todo - how handle deregistration???
				/*$document.bind('click', function() {
				 active = false;
				 })*/
			}
		};

		$document.bind('click', handler);

		scope.$on('$destroy', function() {
			$document.unbind('click', handler);
		})
	}
});