/**
 * @name trail-header
 * @type directive
 *
 * @description Given a trail, constructs the Trail banner
 */
angular.module('clotho.trails').directive('trailHeader', function() {

	var defaultIcon = 'images/trails/trails_splash.jpg';

	return {
		restrict: 'A',
		replace : true,
		templateUrl : 'views/_trails/trailHeader.html',
		scope: {
			trail: '=trailHeader'
		},
		link: function (scope,element,attrs) {
			scope.defaultTrailIcon = defaultIcon;
		}
	}
});