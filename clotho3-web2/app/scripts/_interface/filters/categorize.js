/**
 * todo - major refactor
 * @name categorize
 *
 * @description
 * Organize an array by a field
 *
 * @param [key] {string} the name of the attribute of each object to compare for uniqueness
 if the key is empty, the entire object will be compared
 if the key === false then no filtering will be performed
 *
 * @returns {object} Object with keys matching param used to sort, each containing an array
 */
angular.module('clotho.interface').filter('categorize', function($parse) {
	return function (items, filterOn) {
		if (filterOn === false) {
			return items;
		}

		if ((filterOn || angular.isUndefined(filterOn)) && angular.isArray(items)) {
			var newItems = {},
				get = angular.isString(filterOn) ? $parse(filterOn) : function (item) { return item; };

			var extractValueToCompare = function (item) {
				return angular.isObject(item) ? get(item) : item;
			};

			angular.forEach(items, function (item) {
				var type = extractValueToCompare(item);
				if (!newItems[type])
					newItems[type] = [];
				newItems[type].push(item);
			});

			items = newItems;
		}
		return items;
	}
});