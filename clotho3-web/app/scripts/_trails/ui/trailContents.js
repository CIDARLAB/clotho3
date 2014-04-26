/**
 * @name trail-contents
 * @type directive
 *
 * @description Given a trail, the directive will display its contents (Chapters, Pages)
 */
angular.module('clotho.trails').directive('trailContents', function(Trails) {

	return {
		restrict: 'A',
		replace : true,
		templateUrl : 'views/_trails/trailContents.html',
		scope: {
			trail: '=trailContents'
		},
		link: function (scope,element,attrs) {
			scope.activate = Trails.activate;
			scope.mapIcon = Trails.mapIcon;
		}
	}
});