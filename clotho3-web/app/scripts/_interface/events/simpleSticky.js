/**
 * @ngdoc directive
 * @name simple-sticky
 *
 * @description
 * Simple directive to create a sticky div, which will remain fixed relative to the top of the page.
 *
 * @attr simpleSticky {number} number of pixels from top
 */
angular.module('clotho.interface')
.directive('simpleSticky', function($window) {

	var bodyEl = angular.element($window.document.body);
	var windowEl = angular.element($window);

	return {
		restrict: 'A',
		link: function(scope, element, attrs) {
			var showFromTopSticky = attrs.simpleSticky || '20px',
				showFromTopNormal = element.css('top');

			function parseOffsets () {
				var boundingClientRect = element[0].getBoundingClientRect();
				return boundingClientRect.top + ($window.pageYOffset || $window.document[0].documentElement.scrollTop);
			}

			function checkPosition () {
				//todo - need good way of getting position

				var fromTop = parseOffsets();

				var affixed = true;

				//compare to last time

				if (affixed) {
					element.css({
						top: showFromTopSticky,
						width: "inherit"
					});
				} else {
					element.css({
						top: showFromTopNormal
					});
				}
				element.toggleClass('affix', affixed)
			}


			var checkPositionDebounced = _.debounce(function () {
				checkPosition()
			}, 50);

			//bind listeners

			windowEl.on('scroll', checkPositionDebounced);
			windowEl.on('resize', checkPositionDebounced);

			//remove on destroy

			scope.$on('$destroy', function () {
				windowEl.off('scroll', checkPositionDebounced);
				windowEl.off('resize', checkPositionDebounced);
			});
		}
	}
});