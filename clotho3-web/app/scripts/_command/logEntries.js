/**
 * @name log-entries
 * @example
 *  <div log-entries="myEntries"></div>
 *
 *  where myEntries is an array of entries of the form:
 *  {
 *    class : <string, bootstrap alert CSS class>,
 *    from : <string, 'client' or 'server'>,
 *    text: <string>,
 *    timestamp : Date.now()
 *  }
 */

angular.module('clotho.commandbar')
.directive('logEntries', function() {

	return {
		restrict: 'A',
		templateUrl: 'views/_command/logEntries.html',
		scope: {
			entries: '=logEntries'
		}
	}
});