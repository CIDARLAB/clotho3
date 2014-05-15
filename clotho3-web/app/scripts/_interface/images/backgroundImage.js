angular.module('clotho.interface')
.directive('backgroundImage', function () {
		return {
			restrict: 'A',
			link: function postLink(scope, element, attrs) {
				element.css({
					'background-position' : '50% 50%',
					'background-repeat' : 'none',
					'background-size' : 'cover'
				});

				attrs.$observe('ngSrc', function (newsrc, oldsrc) {
					element.css({
						'background-image': 'url(' + newsrc + ')'
					});
				})
			}
		}
});