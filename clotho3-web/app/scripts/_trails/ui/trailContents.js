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
			trail: '=trailContents',
			current : '=',
			activatePass : '&?activate'
		},
		link: function (scope,element,attrs) {
			scope.activate = function (pos) {
				( angular.isDefined(attrs.activate) && angular.isFunction(scope.activatePass) ) ?
					scope.activatePass({ $position : pos }) :
					Trails.activate.apply(null, arguments);
			};

			scope.mapIcon = Trails.mapIcon;
		}
	}
});