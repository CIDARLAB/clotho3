/**
 * @ngdoc directive
 * @name simple-sticky
 *
 * @description
 * Simple directive to create a sticky div, which will remain fixed relative to the top of the page.
 *
 * Not super efficient (is throttled, but still) so don't use too many on a page
 *
 * @attr simpleSticky {number} number of pixels from top
 *
 * @example
 *
 * <div class="myStickyDiv" simple-sticky="20"></div>
 */
angular.module('clotho.interface')
.directive('simpleSticky', function($window) {

	var windowEl = angular.element($window);

	return {
		restrict: 'A',
		link: function(scope, element, attrs) {
			//70 because of sticky nav at top
			var showFromTopSticky = parseInt(attrs.simpleSticky, 10) || 70,
				showFromTopNormal = element.css('top'),
				startFromTop;

			function setOffset () {
				if (element.hasClass('affix')) {
					//do nothing
				} else {
					var boundingClientRect = element[0].getBoundingClientRect();
					startFromTop = boundingClientRect.top + (angular.isDefined($window.pageYOffset) ? $window.pageYOffset : $window.document[0].documentElement.scrollTop);
				}
			}

			//init
			setOffset();

			function checkPosition () {
				var affixed = ($window.pageYOffset + showFromTopSticky) > startFromTop;

				if (affixed) {
					element.css({
						top: showFromTopSticky + 'px',
						width: "inherit"
					});
				} else {
					element.css({
						top: showFromTopNormal
					});
				}
				element.toggleClass('affix', affixed)
			}


			var checkPositionDebounced = _.throttle(function () {
				checkPosition();
				setOffset();
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