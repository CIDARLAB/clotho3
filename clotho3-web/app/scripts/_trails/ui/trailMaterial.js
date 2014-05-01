/**
 * @name trail-material
 * @type directive
 *
 * @description Given a material object, display it and link etc.
 */
angular.module('clotho.trails').directive('trailMaterial', function(Trails) {
	return {
		restrict: 'A',
		replace : true,
		templateUrl : 'views/_trails/trailMaterial.html',
		scope: {
			material: '=trailMaterial'
		},
		link: function (scope,element,attrs) {
			scope.mapIcon = Trails.mapIcon;
		}
	}
});