angular.module('clotho.interface')
/**
 * @ngdoc directive
 * @name backgroundImage
 *
 * @attr background-contain {Boolean} If truthy, set background-size to contain. Default is cover. (Can only be set initially - no $watch in place)
 *
 * @description
 * Prefer to use ngSrc as way of giving image URL
 *
 * use backgroundImage directly when have base64 etc. which is not whitelisted
 */
.directive('backgroundImage', function () {
		return {
			restrict: 'A',
			link: function postLink(scope, element, attrs) {
				element.css({
					'background-position' : '50% 50%',
					'background-repeat' : 'none',
					'background-size' : (scope.$eval(attrs.backgroundContain) ? 'contain' : 'cover')
				});

				attrs.$observe('ngSrc', function (newsrc, oldsrc) {
					element.css({
						'background-image': 'url(' + newsrc + ')'
					});
				});

				//prefer ngSrc, because of SCE, but doesn't allow for base64
				scope.$watch(function () {
					return attrs.backgroundImage
				}, function (newsrc, oldsrc) {
					element.css({
						'background-image': 'url(' + newsrc + ')'
					});
				});
			}
		}
});