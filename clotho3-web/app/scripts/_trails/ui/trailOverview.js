/**
 * @name trail-overview
 * @type directive
 *
 * @description Given a trail, constructs the Trail overview page
 */
angular.module('clotho.trails')
	.directive('trailOverview', function($sce) {

		return {
			restrict: 'A',
			replace : true,
			templateUrl : 'views/_trails/trailOverview.html',
			scope: {
				trail: '=trailOverview'
			}
		}
	});