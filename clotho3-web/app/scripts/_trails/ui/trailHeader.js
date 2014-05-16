/**
 * @name trail-header
 * @type directive
 *
 * @description Given a trail, constructs the Trail banner
 */
angular.module('clotho.trails')
.directive('trailHeader', function($sce) {

	var defaultIcon = 'images/trails/trails_logo.png';

	return {
		restrict: 'A',
		replace : true,
		templateUrl : 'views/_trails/trailHeader.html',
		scope: {
			trail: '=trailHeader'
		},
		link: function (scope,element,attrs) {
			scope.$watch('trail.icon', function (newsrc) {
				scope.trailIconTrusted = newsrc ? $sce.trustAsResourceUrl(newsrc) : defaultIcon;
			});
		}
	}
});