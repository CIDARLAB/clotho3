/**
 * @name shuffle
 *
 * @example {{myArray | shuffle}}
 *
 * @param [key] {string} if the key === false then no filtering will be performed
 *
 * @return {array}
 *
 */
angular.module('clotho.interface').filter('shuffle',  function() {
	return function (items, filterOn) {
		if (filterOn === false) { return items; }

		if ((filterOn || angular.isUndefined(filterOn)) && angular.isArray(items)) {
			var o = items.slice(0, items.length); // copy
			for(var j, x, i = o.length; i; j = parseInt(Math.random() * i), x = o[--i], o[i] = o[j], o[j] = x);
			items = o;
		}
		return items;
	}
});
