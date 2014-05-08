/**
 * @name clickOutside
 *
 * @description
 * Registers a document click listener which is triggered when clicking on an element not containing this one.
 *
 * Requires the attribute click-outside-active to avoid bloat on each click.
 *
 * Pass in the function directly, without parenthesis.
 *
 * @example
 * <div id="externalClickListener"
        click-outside="myCallback"
        click-outside-active="myBoolean">
 */

angular.module('clotho.interface')
.directive('clickOutside', function($document, $parse) {
	return function(scope, element, attr) {

		//click action will be parsed against the scope on each run, so don't need to watch it
		var clickAction = $parse(attr.clickOutside),
			active;

		//watch to see if active
		scope.$watch(function () {
			return $parse(attr.clickOutsideActive)(scope);
		}, function (newval, oldval) {
			active = !!newval;
		});

		var handler = function (event) {
			if (active) {
				if (!element[0].contains(event.target)) {
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