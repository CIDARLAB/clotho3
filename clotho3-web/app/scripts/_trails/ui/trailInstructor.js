/**
 * @name trail-instructor
 * @type directive
 *
 * @description Given a person object, show them as a course instructor
 */
angular.module('clotho.trails').directive('trailInstructor', function(Clotho) {

	var defaultIcon = 'images/assets/user-default.png';

	return {
		restrict: 'A',
		replace : true,
		templateUrl : 'views/_trails/trailInstructor.html',
		scope: {
			instructorId: '=trailInstructor'
		},
		link: function (scope,element,attrs) {

			scope.instructorIcon = defaultIcon;

			scope.$watch('instructorId', function (newval) {
				Clotho.get(newval).then(function (result) {
					scope.instructor = result;

					scope.instructorIcon = result.icon || defaultIcon;
				});
			});
		}
	}
});